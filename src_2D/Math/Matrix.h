/* Matrix.h:  Header file defining a variety of matrix classes. */

#ifndef _MATRIX_INCLUDED
#define _MATRIX_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cassert>

using namespace std;

/* Include the math and the MV++ vector and matrix macro header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MV_VECTOR_ALL_H_
#include "mvv.h" 
#endif // _MV_VECTOR_ALL_H_

#ifndef _MV_MATRIX_ALL_H_
#include "mvm.h"
#endif // _MV_MATRIX_ALL_H_

/* Define the row vector class, 
   a derived class based on MV++ MV_Vector_double class. */

class ColumnVector;

/********************************************************
 * Class: RowVector                                     *
 *                                                      *
 * Member functions                                     *
 *      size       -- Return vector length.             *
 *      null       -- Return 1 if zero-length vector    *
 *                       and 0 if size > 0.             *
 *      ref        -- Return 1 if vector is a view of   *
 *                             another vector           *
 *                       and 0 otherwise.               *
 *      zero       -- Assign zero vector.               *
 *      one        -- Assign one to every element of    *
 *                    the vector.                       *
 *      unitvec    -- Assign unit vector.               *
 *      abs        -- Absolute value (length) of        *
 *                    vector.                           *
 *      sqr        -- Square of vector.                 *
 *      transpose  -- Return the transpose of vector.   *
 *      norm       -- Return norm of vector.            *
 *      norminf    -- Return infinity norm of vector.   *    
 *                                                      *
 * Member operators                                     *
 *      RV -- a row vector                              *
 *      a  -- a scalar (double)                         *
 *                                                      *
 * RV = RV;                                             *
 * a = RV(i);                                           *
 * a = RV[i];                                           *
 * RV = RV + RV;                                        *
 * RV = RV - RV;                                        *
 * RV = a * RV;                                         *
 * RV = RV * a;                                         *
 * RV = RV / a;                                         *
 * a = RV * RV; (inner product)                         *
 * a = RV * CV; (inner product)                         *
 * RV = +RV;                                            *
 * RV = -RV;                                            *
 * RV += RV;                                            *
 * RV -= RV;                                            *
 * RV == RV;                                            *
 * RV != RV;                                            *
 * cout << RV; (output function)                        *
 * cin  >> RV; (input function)                         *
 *                                                      *
 ********************************************************/
class RowVector: public MV_Vector_double{
  private:
  public:
    /* Creation, copy, and assignment constructors. */
    RowVector(void) { }
    RowVector(unsigned int i)                                    : MV_Vector_double(i)       { /* ... */ }
    RowVector(unsigned int i, const double &xval)                : MV_Vector_double(i, xval) { /* ... */ }
    RowVector(double *xval, unsigned int i)                      : MV_Vector_double(xval, i) { /* ... */ }
    RowVector(const double *xval, unsigned int i)                : MV_Vector_double(xval, i) { /* ... */ }
    RowVector(const MV_Vector_double &V)                         : MV_Vector_double(V)       { /* ... */ }
    RowVector(double *d, unsigned int N, MV_Vector_::ref_type i) : MV_Vector_double(d, N, i) { /* ... */ }
    RowVector(const MV_Vector_double &V, MV_Vector_::ref_type i) : MV_Vector_double(V, i)    { /* ... */ }

    RowVector(const RowVector &RVec, MV_Vector_::ref_type i)     : MV_Vector_double(RVec, i) { /* ... */ }

    /* Destructor. */
    // ~RowVector(void);
    // Use automatically generated destructor.

    /* Assign the zero vector. */
    void zero(void);

    /* Assign the value of one to each element of vector. */
    void one(void);

    /* Assign the unit vector with value of one for Nth element. */
    void unitvec(const int N);

    /* Absolute value (magnitude) of row vector. */
    double abs(void);
    double abs(void) const;
    friend double abs(const RowVector &RVec);

    /* Square of vector. */
    double sqr(void);
    double sqr(void) const;
    friend double sqr(const RowVector &RVec);

    /* Return the transpose of the row vector. */
    ColumnVector transpose(void);
    ColumnVector transpose(void) const;
    friend ColumnVector transpose(const RowVector &RVec);

    /* Various norms of row vector. */
    double norm(void);
    double norm(void) const;
    friend double norm(const RowVector &RVec);

    double norm(const int p);
    double norm(const int p) const;
    friend double norm(const RowVector &RVec, const int p);

    double norminf(void);
    double norminf(void) const;
    friend double norminf(const RowVector &RVec);

    /* Assignment operator. */
    // RowVector operator = (const RowVector &RVec);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend RowVector operator /(const RowVector &RVec1, const double &a);
    friend double operator *(const RowVector &RVec1, const RowVector &RVec2);
    friend double operator *(const RowVector &RVec, const ColumnVector &CVec);

    /* Relational operators. */
    friend int operator ==(const RowVector &RVec1, const RowVector &RVec2);
    friend int operator !=(const RowVector &RVec1, const RowVector &RVec2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const RowVector &RVec);
    friend istream &operator >> (istream &in_file, RowVector &RVec);

};

/********************************************************
 * RowVector::zero -- Assign zero vector.               *
 ********************************************************/
inline void RowVector::zero(void) {
   int i, m; m = size(); for ( i = 0; i <= m-1; ++i ) p_[i]=ZERO;
}

/********************************************************
 * RowVector::one -- Assign one to each element.        *
 ********************************************************/
inline void RowVector::one(void) {
   int i, m; m = size(); for ( i = 0; i <= m-1; ++i ) p_[i]=ONE;
}

/********************************************************
 * RowVector::unitvec -- Assign unit vector.            *
 ********************************************************/
inline void RowVector::unitvec(const int N) {
   int i, m; m = size(); assert( N >= 0 && N <= m-1 );
   for ( i = 0; i <= m-1; ++i ) p_[i]=ZERO; p_[N]=ONE;
}

/********************************************************
 * RowVector::abs -- Absolute value of a vector.        *
 ********************************************************/
inline double RowVector::abs(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (sqrt(xx));
}

inline double RowVector::abs(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (sqrt(xx));
}

inline double abs(const RowVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += RVec(i)*RVec(i);
   return (sqrt(xx));
}

/********************************************************
 * RowVector::sqr -- Square of a vector.                *
 ********************************************************/
inline double RowVector::sqr(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (xx);
}

inline double RowVector::sqr(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (xx);
}

inline double sqr(const RowVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += RVec(i)*RVec(i);
   return (xx);
}

/********************************************************
 * RowVector::norm -- Norm of a vector.                 *
 ********************************************************/
inline double RowVector::norm(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(p_[i]);
   return (xx);
}

inline double RowVector::norm(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(p_[i]);
   return (xx);
}

inline double norm(const RowVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(RVec(i));
   return (xx);
}

inline double RowVector::norm(const int p) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(p_[i]), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double RowVector::norm(const int p) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(p_[i]), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double norm(const RowVector &RVec, const int p) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(RVec(i)), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double RowVector::norminf(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(p_[i]));
   return (xx);
}

inline double RowVector::norminf(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(p_[i]));
   return (xx);
}

inline double norminf(const RowVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(RVec(i)));
   return (xx);
}

/********************************************************
 * RowVector -- Binary arithmetic operators.            *
 ********************************************************/
inline RowVector operator /(const RowVector &RVec, const double &a) {
   int i, m; m = RVec.size(); double xx; xx = ONE/a; RowVector rv(m);
   for ( i = 0; i <= m-1; ++i ) {
      rv(i) = xx*RVec(i);
   } /* endfor */
   return (rv);
}

inline double operator *(const RowVector &RVec1, const RowVector &RVec2) {
   return (dot(RVec1, RVec2));
}

/********************************************************
 * RowVector -- Relational operators.                   *
 ********************************************************/
inline int operator ==(const RowVector &RVec1, const RowVector &RVec2) {
   if (RVec1.null() || RVec2.null()) return(0);
   int i, equal; i = 0; equal=1;
   if (RVec1.size() != RVec2.size()) equal = 0;
   while (equal) {
      if (RVec1(i) != RVec2(i)) equal = 0;
      if (i == (int)RVec1.size()-1) break;
      i = i + 1;
   } /* endwhile */
   return (equal);
}

inline int operator !=(const RowVector &RVec1, const RowVector &RVec2) {
   if (RVec1.null() || RVec2.null()) return(1);
   int i, not_equal; i = 0; not_equal=1;
   if (RVec1.size() == RVec2.size()) {
      while (not_equal) {
         if (RVec1(i) == RVec2(i)) not_equal = 0;
         if (i ==  (int)RVec1.size()-1) break;
         i = i + 1;
      } /* endwhile */
   } /* endif */
   return (not_equal);
}

/********************************************************
 * RowVector -- Input-output operators.                 *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const RowVector &RVec) {
   int i, m; m = RVec.size();
   out_file.setf(ios::scientific);
   for ( i = 0; i <= m-1; ++i ) {
       out_file << " " << RVec(i);
   } /* endfor */
   out_file.unsetf(ios::scientific); out_file << "\n";
   return (out_file);
}

inline istream &operator >> (istream &in_file, RowVector &RVec) {
   int i, m; m = RVec.size();
   in_file.setf(ios::skipws);
   for ( i = 0; i <= m-1; ++i ) {
       in_file >> RVec(i);
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/* Define the column vector class, 
   a derived class based on MV++ MV_Vector_double class. */

/********************************************************
 * Class: ColumnVector                                  *
 *                                                      *
 * Member functions                                     *
 *      size       -- Return vector length.             *
 *      null       -- Return 1 if zero-length vector    *
 *                       and 0 if size > 0.             *
 *      ref        -- Return 1 if vector is a view of   *
 *                             another vector           *
 *                       and 0 otherwise.               *
 *      zero       -- Assign zero vector.               *
 *      one        -- Assign one to every element of    *
 *                    the vector.                       *
 *      unitvec    -- Assign unit vector.               *
 *      abs        -- Absolute value (length) of        *
 *                    vector.                           *
 *      sqr        -- Square of vector.                 *
 *      transpose  -- Return the transpose of vector.   *
 *      norm       -- Return norm of vector.            *
 *      norminf    -- Return infinity norm of vector.   *
 *                                                      *
 * Member operators                                     *
 *      CV -- a column vector                           *
 *      RV -- a row vector                              *
 *      a  -- a scalar (double)                         *
 *                                                      *
 * CV = CV;                                             *
 * a = CV(i);                                           *
 * a = CV[i];                                           *
 * CV = CV + CV;                                        *
 * CV = CV - CV;                                        *
 * CV = a * CV;                                         *
 * CV = CV * a;                                         *
 * CV = CV / a;                                         *
 * a = CV * CV; (inner product)                         *
 * a = RV * CV; (inner product)                         *
 * CV = +CV;                                            *
 * CV = -CV;                                            *
 * CV += CV;                                            *
 * CV -= CV;                                            *
 * CV == CV;                                            *
 * CV != CV;                                            *
 * cout << CV; (output function)                        *
 * cin  >> CV; (input function)                         *
 *                                                      *
 ********************************************************/
class ColumnVector: public MV_Vector_double{
  private:
  public:
    /* Creation, copy, and assignment constructors. */
    ColumnVector(void)                                              : MV_Vector_double() { }
    ColumnVector(unsigned int i)                                    : MV_Vector_double(i)       { /* ... */ }
    ColumnVector(unsigned int i, const double &xval)                : MV_Vector_double(i, xval) { /* ... */ }
    ColumnVector(double *xval, unsigned int i)                      : MV_Vector_double(xval, i) { /* ... */ }
    ColumnVector(const double *xval, unsigned int i)                : MV_Vector_double(xval, i) { /* ... */ }
    ColumnVector(const MV_Vector_double &V)                         : MV_Vector_double(V)       { /* ... */ }
    ColumnVector(double *d, unsigned int N, MV_Vector_::ref_type i) : MV_Vector_double(d, N, i) { /* ... */ }
    ColumnVector(const MV_Vector_double &V, MV_Vector_::ref_type i) : MV_Vector_double(V, i)    { /* ... */ }

    ColumnVector(const ColumnVector &RVec, MV_Vector_::ref_type i)  : MV_Vector_double(RVec, i) { /* ... */ }

    /* Destructor. */
    // ~ColumnVector(void);
    // Use automatically generated destructor.

    /* Assign the zero vector. */
    void zero(void);

    /* Assign the value of one to each element of vector. */
    void one(void);

    /* Assign the unit vector with value of one for Nth element. */
    void unitvec(const int N);

    /* Absolute value (magnitude) of row vector. */
    double abs(void);
    double abs(void) const;
    friend double abs(const ColumnVector &CVec);

    /* Square of vector. */
    double sqr(void);
    double sqr(void) const;
    friend double sqr(const ColumnVector &CVec);

    /* Return the transpose of the column vector. */
    RowVector transpose(void);
    RowVector transpose(void) const;
    friend RowVector transpose(const ColumnVector &CVec);
    
    /* Various norms of column vector. */
    double norm(void);
    double norm(void) const;
    friend double norm(const ColumnVector &CVec);

    double norm(const int p);
    double norm(const int p) const;
    friend double norm(const ColumnVector &CVec, const int p);

    double norminf(void);
    double norminf(void) const;
    friend double norminf(const ColumnVector &CVec);

    /* Assignment operator. */
    // ColumnVector operator = (const ColumnVector &CVec);
    // Use automatically generated assignment operator.

    /* Binary arithmetic operators. */
    friend ColumnVector operator /(const ColumnVector &CVec1, const double &a);
    friend double operator *(const ColumnVector &CVec1, const ColumnVector &CVec2);
    friend double operator *(const RowVector &RVec, const ColumnVector &CVec);
    friend ColumnVector operator *(const double &a, const ColumnVector &CVec1);
    /* term-wise product */
    ColumnVector  operator ^(const ColumnVector &CVec2);

    /* Relational operators. */
    friend int operator ==(const ColumnVector &CVec1, const ColumnVector &CVec2);
    friend int operator !=(const ColumnVector &CVec1, const ColumnVector &CVec2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const ColumnVector &CVec1);
    friend istream &operator >> (istream &in_file, ColumnVector &CVec1);

};

/********************************************************
 * ColumnVector::zero -- Assign zero vector.            *
 ********************************************************/
inline void ColumnVector::zero(void) {
   int i, m; m = size(); for ( i = 0; i <= m-1; ++i ) p_[i]=ZERO;
}

/********************************************************
 * ColumnVector::one -- Assign one to each element.     *
 ********************************************************/
inline void ColumnVector::one(void) {
   int i, m; m = size(); for ( i = 0; i <= m-1; ++i ) p_[i]=ONE;
}

/********************************************************
 * ColumnVector::unitvec -- Assign unit vector.         *
 ********************************************************/
inline void ColumnVector::unitvec(const int N) {
   int i, m; m = size(); assert( N >= 0 && N <= m-1 );
   for ( i = 0; i <= m-1; ++i ) p_[i]=ZERO; p_[N]=ONE;
}

/********************************************************
 * ColumnVector::abs -- Absolute value of a vector.     *
 ********************************************************/
inline double ColumnVector::abs(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (sqrt(xx));
}

inline double ColumnVector::abs(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (sqrt(xx));
}

inline double abs(const ColumnVector &CVec) {
   int i, m; m = CVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += CVec(i)*CVec(i);
   return (sqrt(xx));
}

/********************************************************
 * ColumnVector::sqr -- Square of a vector.             *
 ********************************************************/
inline double ColumnVector::sqr(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (xx);
}

inline double ColumnVector::sqr(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += p_[i]*p_[i];
   return (xx);
}

inline double sqr(const ColumnVector &CVec) {
   int i, m; m = CVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += CVec(i)*CVec(i);
   return (xx);
}

/********************************************************
 * ColumnVector::norm -- Norm of a vector.              *
 ********************************************************/
inline double ColumnVector::norm(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(p_[i]);
   return (xx);
}

inline double ColumnVector::norm(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(p_[i]);
   return (xx);
}

inline double norm(const ColumnVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += fabs(RVec(i));
   return (xx);
}

inline double ColumnVector::norm(const int p) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(p_[i]), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double ColumnVector::norm(const int p) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(p_[i]), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double norm(const ColumnVector &RVec, const int p) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx += pow(fabs(RVec(i)), double(p));
   return (pow(xx, ONE/double(p)));
}

inline double ColumnVector::norminf(void) {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(p_[i]));
   return (xx);
}

inline double ColumnVector::norminf(void) const {
   int i, m; m = size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(p_[i]));
   return (xx);
}

inline double norminf(const ColumnVector &RVec) {
   int i, m; m = RVec.size(); double xx; xx= ZERO;
   for ( i = 0; i <= m-1; ++i ) xx = max(xx, fabs(RVec(i)));
   return (xx);
}

/********************************************************
 * ColumnVector -- Binary arithmetic operators.         *
 ********************************************************/
inline ColumnVector operator /(const ColumnVector &CVec, const double &a) {
   int i, m; m = CVec.size(); double xx; xx = ONE/a; ColumnVector cv(m);
   for ( i = 0; i <= m-1; ++i ) {
      cv(i) = xx*CVec(i);
   } /* endfor */
   return (cv);
}

inline ColumnVector ColumnVector::operator ^(const ColumnVector &CVec2) {
  int i, m; m = size(); ColumnVector cv(m);
  for ( i = 0; i <= m-1; ++i ) {
    cv(i) = (*this)(i)*CVec2(i);
  } /* endfor */
  return (cv);
}

inline ColumnVector operator *(const double &a, const ColumnVector &CVec) {
   int i, m; m = CVec.size(); ColumnVector cv(m);
   for ( i = 0; i <= m-1; ++i ) {
      cv(i) = a*CVec(i);
   } /* endfor */
   return (cv);
}

inline double operator *(const ColumnVector &CVec1, const ColumnVector &CVec2) {
   return (dot(CVec1, CVec2));
}

/********************************************************
 * ColumnVector -- Relational operators.                *
 ********************************************************/
inline int operator ==(const ColumnVector &CVec1, const ColumnVector &CVec2) {
   if (CVec1.null() || CVec2.null()) return(0);
   int i, equal; i = 0; equal=1;
   if (CVec1.size() != CVec2.size()) equal = 0;
   while (equal) {
      if (CVec1(i) != CVec2(i)) equal = 0;
      if (i ==  (int)CVec1.size()-1) break;
      i = i + 1;
   } /* endwhile */
   return (equal);
}

inline int operator !=(const ColumnVector &CVec1, const ColumnVector &CVec2) {
   if (CVec1.null() || CVec2.null()) return(1);
   int i, not_equal; i = 0; not_equal=1;
   if (CVec1.size() == CVec2.size()) {
      while (not_equal) {
         if (CVec1(i) == CVec2(i)) not_equal = 0;
         if (i ==  (int)CVec1.size()-1) break;
         i = i + 1;
      } /* endwhile */
   } /* endif */
   return (not_equal);
}

/********************************************************
 * ColumnVector -- Input-output operators.              *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const ColumnVector &CVec) {
  int i, m; m = CVec.size();
  out_file.setf(ios::scientific);
  for ( i = 0; i <= m-1; ++i ) {
      out_file << " " << CVec(i) << "\n";
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, ColumnVector &CVec) {
  int i, m; m = CVec.size();
  in_file.setf(ios::skipws);
  for ( i = 0; i <= m-1; ++i ) {
      in_file >> CVec(i);
  } /* endfor */
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * RowVector::transpose -- Return vector transpose.     *
 ********************************************************/
inline ColumnVector RowVector::transpose(void) {
   return (ColumnVector(p_, dim_));
}

inline ColumnVector RowVector::transpose(void) const {
   return (ColumnVector(p_, dim_));
}

inline ColumnVector transpose(const RowVector &RVec) {
   return (ColumnVector(RVec));
}

/********************************************************
 * RowVector -- Binary arithmetic operators.            *
 ********************************************************/
inline double operator *(const RowVector &RVec, const ColumnVector &CVec) {
   return (dot(RVec, CVec));
}

/********************************************************
 * ColumnVector::transpose -- Return vector transpose.  *
 ********************************************************/
inline RowVector ColumnVector::transpose(void) {
   return (RowVector(p_, dim_));
}

inline RowVector ColumnVector::transpose(void) const {
   return (RowVector(p_, dim_));
}

inline RowVector transpose(const ColumnVector &CVec) {
   return (RowVector(CVec));
}

/* Define the regular dense matrix class, 
   a derived class based on MV++ MV_ColMat_double class. */

/********************************************************
 * Class: DenseMatrix                                   *
 *                                                      *
 * Member functions                                     *
 *      size       -- Return matrix dimensions (mxn).   *
 *      lda        -- Return leading dimension of       *
 *                    matrix.                           *
 *      ref        -- Return 1 if matrix is a view of   *
 *                             another C/C++ array      *
 *                       and 0 otherwise.               *
 *      null       -- Return 1 if zero-length matrix    *
 *                       and 0 if size > 0.             *
 *      zero       -- Assign zero matrix.               *
 *      identity   -- Assign identity matrix.           *
 *      trace      -- Return the trace of matrix.       *
 *      transpose  -- Return the transpose of matrix.   *
 *      permute_col-- Return the matrix having the      *
 *                    col_1 permuted with col_2         *
 *      permute_row-- Return the matrix having the      *
 *                    row_1 permuted with row_2         *
 * pseudo_inverse_override -- Computes the              *
 *                        pseudo-inverse of             *
 *                        a dense matrix MxN and writes *
 *                        it on top of the initial      *
 *                        matrix.                       *
 *                                                      *
 * pseudo_inverse -- Returns the pseudo-inverse of      *
 *                   a dense matrix MxN.                *
 *                                                      *
 *  eigenvalues   -- Return a vector containnig the     *
 *                   eigenvalues of an NxN matrix.      *
 *                                                      *
 *    inverse     -- Return the inverse of a matrix     *
 *                   using LAPACK's dgetri function.    *
 * Member operators                                     *
 *      M  -- a regular mxn dense matrix                *
 *      RV -- a row vector                              *
 *      CV -- a column vector                           *
 *      a  -- a scalar (double)                         *
 *                                                      *
 * M = M;                                               *
 * a = M(i,j);                                          *
 * RV = M[i];                                           *
 * a = M[i][j];                                         * 
 * M = M + M;                                           *
 * M = M - M;                                           *
 * M = M * M;                                           *
 * CV = M * CV;                                         *
 * M = CV*RV;                                           *
 * M = a * M;                                           *
 * M = M * a;                                           *
 * M = M / a;                                           *
 * M = +M;                                              *
 * M = -M;                                              *
 * M += M;                                              *
 * M -= M;                                              *
 * M == M;                                              *
 * M != M;                                              *
 * cout << M; (output function)                         *
 * cin  >> M; (input function)                          *
 *                                                      *
 ********************************************************/
class DenseMatrix: public MV_ColMat_double{
  private:
    static RowVector temp_RVec;
  public:
    /* Creation, copy, and assignment constructors. */
    DenseMatrix(void) : MV_ColMat_double() { }
    DenseMatrix(unsigned int m, unsigned int n) : 
       MV_ColMat_double(m, n) { /* ... */ } 
    DenseMatrix(unsigned int m, unsigned int n, const double &x) : 
       MV_ColMat_double(m, n, x) { /* ... */ }   
    DenseMatrix(double *x, unsigned int m, unsigned int n) : 
       MV_ColMat_double(x, m, n) { /* ... */ }
    DenseMatrix(double *x, unsigned int m, unsigned int n, unsigned int lda) : 
       MV_ColMat_double(x, m, n, lda) { /* ... */ }
    DenseMatrix(MV_ColMat_double &A, MV_Matrix_::ref_type i) : 
       MV_ColMat_double(A, i) { /* ... */ }
    DenseMatrix(double *x, unsigned int m, unsigned int n, MV_Matrix_::ref_type i) : 
       MV_ColMat_double(x, m, n, i) { /* ... */ }
    DenseMatrix(double *x, unsigned int m, unsigned int n, unsigned int lda, MV_Matrix_::ref_type i) : 
       MV_ColMat_double(x, m, n, lda, i) { /* ... */ }
    DenseMatrix(const MV_ColMat_double &A) :
       MV_ColMat_double(A) { /* ... */ }

    DenseMatrix(DenseMatrix &A, MV_Matrix_::ref_type i) : MV_ColMat_double(A, i) { /* ... */ }

    /* Destructor. */
    // ~DenseMatrix(void);
    // Use automatically generated destructor.

    /* Find out if the matrix has no elements. */
    inline int null(void) const {return dim0_== 0 && dim1_==0;}

    /* Obtain the dimensions of the mxn matrix. */
    int get_n(void) {return dim0_;}
    int get_m(void) {return dim1_;}

    /* Assign the zero matrix. */
    void zero(void);

    /* Assign the identity matrix. */
    void identity(void);

    /* Compute the trace of the matrix. */
    double trace(void);
    double trace(void) const;
    friend double trace(const DenseMatrix &M);

    /* Return the transpose of the matrix. */
    DenseMatrix transpose(void);
    DenseMatrix transpose(void) const;
    friend DenseMatrix transpose(const DenseMatrix &M);
    
    /* Get the pseudo-inverse for an over/under-determined matrix */
    DenseMatrix pseudo_inverse(void) const;
    friend DenseMatrix pseudo_inverse(const DenseMatrix &A);
    void pseudo_inverse_override(void);
    friend void pseudo_inverse_override(DenseMatrix &A);

    /* Get a vector containing the eigenvalues of an NxN matrix */
    ColumnVector eigenvalues(void) const;
    ColumnVector eigenvalues_overwrite(void);

    /* Get matrix inverse */
    DenseMatrix inverse(void) const;
    void inverse_overwrite(void);

    /* Compute the Frobenious norm of the matrix ||A|| = sqrt(SUM_i SUM_j (A(i,j)^2) ) */
    double NormFro (void) const;

    /* Return the matrix having row1_ permuted with row2_ */
    void permute_row(const int col1_, const int col2_);
    DenseMatrix permute_row(const int row1_, const int row2_) const;
    friend DenseMatrix permute_row(const DenseMatrix &M, const int row1_,
				   const int row2_);

    /* Return the matrix having col1_ permuted with col2_ */
    void permute_col(const int col1_, const int col2_);
    DenseMatrix permute_col(const int col1_, const int col2_) const;
    friend DenseMatrix permute_col(const DenseMatrix &M, const int col1_,
				   const int col2_);
    
    /*! @brief Make matrix M part of the current matrix starting from (StartRow,StartCol) element */
    void incorporate_matrix(const int & StartRow, const int & StartCol, 
			    const DenseMatrix &M);

    /* Assignment operator. */
    // DenseMatrix operator = (const DenseMatrix &M1);
    // Use automatically generated assignment operator.

    /* Index operator. */
    RowVector &operator[](int index) {
       assert( index >= 0 && index <= dim0_-1); int j; 
       temp_RVec = RowVector(dim1_, ZERO);
       for ( j = 0 ; j <= dim1_-1; ++j ) {
          temp_RVec(j) = v_(j*lda_+index);
       } /* endfor */
       return (temp_RVec);
    }
    
    const RowVector &operator[](int index) const {
       assert( index >= 0 && index <= dim0_-1); int j; 
       temp_RVec = RowVector(dim1_, ZERO);
       for ( j = 0 ; j <= dim1_-1; ++j ) {
          temp_RVec(j) = v_(j*lda_+index);
       } /* endfor */
       return (temp_RVec);
    }

    /* Unary arithmetic operators. */
    friend DenseMatrix operator +(const DenseMatrix &M);
    friend DenseMatrix operator -(const DenseMatrix &M);

    /* Binary arithmetic operators. */
    friend DenseMatrix operator +(const DenseMatrix &M1, const DenseMatrix &M2);
    friend DenseMatrix operator -(const DenseMatrix &M1, const DenseMatrix &M2);
    friend DenseMatrix operator *(const DenseMatrix &M1, const DenseMatrix &M2);
    friend ColumnVector operator *(const DenseMatrix &M, const ColumnVector &CVec);
    friend DenseMatrix operator *(const ColumnVector &CVec, const RowVector &RVec);
    friend DenseMatrix operator *(const DenseMatrix &M, const double &a);
    friend DenseMatrix operator *(const double &a, const DenseMatrix &M);
    friend DenseMatrix operator /(const DenseMatrix &M, const double &a);

    /* Shortcut arithmetic operators. */
    friend DenseMatrix &operator +=(DenseMatrix &M1, const DenseMatrix &M2);
    friend DenseMatrix &operator -=(DenseMatrix &M1, const DenseMatrix &M2);
    friend DenseMatrix &operator *=(DenseMatrix &M1, double c);
    
    /* Relational operators. */
    friend int operator ==(const DenseMatrix &M1, const DenseMatrix &M2);
    friend int operator !=(const DenseMatrix &M1, const DenseMatrix &M2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const DenseMatrix &M);
    friend istream &operator >> (istream &in_file, DenseMatrix &M);

};

/********************************************************
 * DenseMatrix::zero -- Assign zero matrix.             *
 ********************************************************/
inline void DenseMatrix::zero(void) {
   int i, j, k;
   for ( i = 0; i <= dim0_-1; ++i ) {
      for ( j = 0 ; j <= dim1_-1; ++j ) {
         k=j*lda_+i; v_(k)=ZERO;
      } /* endfor */
   } /* endfor */
}

/********************************************************
 * DenseMatrix::identity -- Assign identity matrix.     *
 ********************************************************/
inline void DenseMatrix::identity(void) {
   int i, j, k;
   for ( i = 0; i <= dim0_-1; ++i ) {
      for ( j = 0 ; j <= dim1_-1; ++j ) {
         k=j*lda_+i; v_(k)=double(i==j)*ONE;
      } /* endfor */
   } /* endfor */
}

/********************************************************
 * DenseMatrix::trace -- Trace of a matrix.             *
 ********************************************************/
inline double DenseMatrix::trace(void) {
   int i, k; double xx; xx = ZERO;
   for ( i = 0; i <= min(int(dim0_-1),int(dim1_-1)); ++i ) {
      k=i*(lda_+1); xx+=v_(k);
   } /* endfor */
   return (xx);
}

inline double DenseMatrix::trace(void) const {
   int i, k; double xx; xx = ZERO;
   for ( i = 0; i <= min(int(dim0_-1),int(dim1_-1)); ++i ) {
      k=i*(lda_+1); xx+=v_(k);
   } /* endfor */
   return (xx);
}

inline double trace(const DenseMatrix &M) {
   int i; double xx; xx = ZERO;
   for ( i = 0; i <= min(int(M.dim0_-1),int(M.dim1_-1)); ++i ) {
      xx+=M(i,i);
   } /* endfor */
   return (xx);
}

/********************************************************
 * DenseMatrix::transpose -- Transpose of a matrix.     *
 ********************************************************/

inline DenseMatrix DenseMatrix::transpose(void) {
   int i, j, k; DenseMatrix MT(dim1_,dim0_);
   for ( j = 0 ; j <= dim1_-1; ++j ) {
      for ( i = 0; i <= dim0_-1; ++i ) {
          k=j*lda_+i; MT(j,i)=v_(k);
      } /* endfor */
   } /* endfor */
   return (MT);
}

inline DenseMatrix DenseMatrix::transpose(void) const {
   int i, j, k; DenseMatrix MT(dim1_,dim0_);
   for ( j = 0 ; j <= dim1_-1; ++j ) {
      for ( i = 0; i <= dim0_-1; ++i ) {
          k=j*lda_+i; MT(j,i)=v_(k);
      } /* endfor */
   } /* endfor */
   return (MT);
}

inline DenseMatrix transpose(const DenseMatrix &M) {
   int i, j; DenseMatrix MT(M.dim1_,M.dim0_);
   for ( j = 0 ; j <= M.dim1_-1; ++j ) {
      for ( i = 0; i <= M.dim0_-1; ++i ) {
          MT(j,i) = M(i,j);
      } /* endfor */
   } /* endfor */
   return (MT);
}

/*********************************************************************
 * pseudo_inverse() -- Return the pseudo-inverse of the imput matrix *
 *******************************************************************/
inline DenseMatrix DenseMatrix::pseudo_inverse(void) const{
  DenseMatrix Copy(*this);
  Copy.pseudo_inverse_override(); /* compute pseudo-inverse for the Copy */
  return Copy;
}

/****************************************************************************
 * friend pseudo_inverse() -- Return the pseudo-inverse of the imput matrix *
 ***************************************************************************/
inline DenseMatrix pseudo_inverse(const DenseMatrix &A){
  return A.pseudo_inverse();
}

/*******************************************************************
 * friend pseudo_inverse_override() -- Compute the pseudo-inverse  *
 *                                     and write it on top of the  *
 *                                     original matrix.            *
 *******************************************************************/
inline void pseudo_inverse_override(DenseMatrix &A){
  return A.pseudo_inverse_override();
}

/*************************************************************
 * eigenvalues -- Return a ColumnVector containng the        *
 *                eigenvalues of an NxN matrix.              *
 *************************************************************/
inline ColumnVector DenseMatrix::eigenvalues(void) const{
  DenseMatrix Copy(*this);
  return Copy.eigenvalues_overwrite();
}

/*************************************************************
 * eigenvalues -- Return a ColumnVector containng the        *
 *                eigenvalues of an NxN matrix.              *
 *************************************************************/
inline DenseMatrix DenseMatrix::inverse(void) const {
  DenseMatrix Copy(*this);
  Copy.inverse_overwrite();
  return Copy;
}

/****************************************************************
 * DenseMatrix::NormFro -- Permutation of row1_ with row2_. *
 ***************************************************************/

inline double DenseMatrix::NormFro(void) const{
  double Norm(0.0);
  int i,j,k;

  for ( j = 0; j <=  (int)dim(1)-1; ++j ) {
    for ( i = 0; i <=  (int)dim(0)-1; ++i ) {
      k=j*lda()+i;
      Norm += v_(k)*v_(k);
    } // endfor
  } // endfor
  
  return sqrt(Norm);
}

/****************************************************************
 * DenseMatrix::permute_row -- Permutation of row1_ with row2_. *
 ***************************************************************/
inline void DenseMatrix::permute_row(const int row1_, const int row2_)
{
  double temp;
  for (int j=0; j< (int)dim(1); ++j){
    temp = operator()(row1_,j);
    operator()(row1_,j) = operator()(row2_,j);
    operator()(row2_,j) = temp;
  }
}

inline DenseMatrix DenseMatrix::permute_row(const int row1_, const int row2_) const {
  DenseMatrix MT;

  MT = *this;
  MT.permute_row(row1_,row2_);
  return MT;
}

inline DenseMatrix permute_row(const DenseMatrix &M, const int row1_,
			       const int row2_) {
  return M.permute_row(row1_,row2_);
}

/****************************************************************
 * DenseMatrix::permute_col -- Permutation of col1_ with col2_. *
 ***************************************************************/
inline void DenseMatrix::permute_col(const int col1_, const int col2_)
{
  double temp;
  for (int i=0; i< (int)dim(0); ++i){
    temp = operator()(i,col1_);
    operator()(i,col1_) = operator()(i,col2_);
    operator()(i,col2_) = temp;
  }
}

inline DenseMatrix DenseMatrix::permute_col(const int col1_, const int col2_) const {
  DenseMatrix MT;

  MT = *this;
  MT.permute_col(col1_,col2_);
  return MT;
}

inline DenseMatrix permute_col(const DenseMatrix &M, const int col1_,
			       const int col2_) {

  return M.permute_col(col1_,col2_);
}

/********************************************************
 * DenseMatrix -- Unary arithmetic operators.           *
 ********************************************************/
inline DenseMatrix operator +(const DenseMatrix &M) {
   return (M);
}

inline DenseMatrix operator -(const DenseMatrix &M) {
   DenseMatrix MM(M.dim0_,M.dim1_); MM.v_ = (-ONE)*M.v_;
   return (MM);
}

/********************************************************
 * DenseMatrix -- Binary arithmetic operators.          *
 ********************************************************/
inline DenseMatrix operator +(const DenseMatrix &M1, const DenseMatrix &M2) {
   assert(M1.dim0_ == M2.dim0_ && M1.dim1_ == M2.dim1_);
   DenseMatrix M(M1.dim0_,M1.dim1_); M.v_ = M1.v_ + M2.v_;
   return (M);
}

inline DenseMatrix operator -(const DenseMatrix &M1, const DenseMatrix &M2) {
   assert(M1.dim0_ == M2.dim0_ && M1.dim1_ == M2.dim1_);
   DenseMatrix M(M1.dim0_,M1.dim1_); M.v_ = M1.v_ - M2.v_;
   return (M);
}

inline DenseMatrix operator *(const DenseMatrix &M1, const DenseMatrix &M2) {
  assert(M1.dim1_ == M2.dim0_);
   int i, j, k; DenseMatrix MP(M1.dim0_,M2.dim1_,ZERO);
   for ( i = 0; i <= MP.dim0_-1; ++i ) {
      for ( j = 0 ; j <= MP.dim1_-1; ++j ) {
	  for ( k = 0 ; k <= M1.dim1_-1; ++k ) { MP(i,j) += M1(i,k)*M2(k,j); }
      } /* endfor */
   } /* endfor */
   return (MP);
}

inline ColumnVector operator *(const DenseMatrix &M, const ColumnVector &CVec) {
   assert(M.dim1_ ==  (int)CVec.size());
   int i, j; ColumnVector cv(M.dim0_, ZERO);
   for ( i = 0; i <= M.dim0_-1; ++i ) {
      for ( j = 0 ; j <= M.dim1_-1; ++j ) cv(i)+=M(i,j)*CVec(j);
   } /* endfor */
   return (cv);
}

inline DenseMatrix operator *(const ColumnVector &CVec, const RowVector &RVec){
   int i, j; DenseMatrix M(CVec.size(), RVec.size());
   for ( i = 0; i <= M.dim0_-1; i++ ) {
      for ( j = 0 ; j <= M.dim1_-1; j++ )
	M(i,j) = CVec(i)*RVec(j);
   } /* endfor */
   return (M);
}


inline DenseMatrix operator *(const DenseMatrix &M, const double &a) {
   DenseMatrix Ma(M.dim0_,M.dim1_); Ma.v_ = a*M.v_;
   return (Ma);
}

inline DenseMatrix operator *(const double &a, const DenseMatrix &M) {
   DenseMatrix Ma(M.dim0_,M.dim1_); Ma.v_ = a*M.v_;
   return (Ma);
}

inline DenseMatrix operator /(const DenseMatrix &M, const double &a) {
   double xx; xx = ONE/a; DenseMatrix Ma(M.dim0_,M.dim1_); Ma.v_ = xx*M.v_;
   return (Ma);
}

/********************************************************
 * DenseMatrix -- Shortcut arithmetic operators.        *
 ********************************************************/
inline DenseMatrix &operator +=(DenseMatrix &M1, const DenseMatrix &M2) {
   assert(M1.dim0_ == M2.dim0_ && M1.dim1_ == M2.dim1_);
   M1.v_ += M2.v_;
   return (M1);
}

inline DenseMatrix &operator -=(DenseMatrix &M1, const DenseMatrix &M2) {
   assert(M1.dim0_ == M2.dim0_ && M1.dim1_ == M2.dim1_);
   M1.v_ -= M2.v_;
   return (M1);
}

inline DenseMatrix &operator *=(DenseMatrix &M1, double c) {
   M1.v_ *= c;
   return (M1);
}

/********************************************************
 * DenseMatrix -- Relational operators.                 *
 ********************************************************/
inline int operator ==(const DenseMatrix &M1, const DenseMatrix &M2) {
   if (M1.dim0_ != M2.dim0_ || M1.dim1_ != M2.dim1_) return(0);
   int i , j, equal; equal = 1;
   for ( i = 0; i <= M1.dim0_-1; ++i ) {
      for ( j = 0 ; j <= M1.dim1_-1; ++j ) {
	 if (M1(i,j) != M2(i,j)) {
            equal = 0;
            break;
         } /* endif */
      } /* endfor */
      if (!equal) break;
   } /* endfor */
   return (equal);
}

inline int operator !=(const DenseMatrix &M1, const DenseMatrix &M2) {
   if (M1.dim0_ != M2.dim0_ || M1.dim1_ != M2.dim1_) return(1);
   int i , j, not_equal; not_equal = 1;
   for ( i = 0; i <= M1.dim0_-1; ++i ) {
      for ( j = 0 ; j <= M1.dim1_-1; ++j ) {
	 if (M1(i,j) == M2(i,j)) {
            not_equal = 0;
            break;
         } /* endif */
      } /* endfor */
      if (!not_equal) break;
   } /* endfor */
   return (not_equal);
}

/********************************************************
 * DenseMatrix -- Input-output operators.               *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const DenseMatrix &M) {
  int i;
  out_file.setf(ios::scientific);
  for ( i = 0 ; i <= M.dim0_-1; ++i ) {
      out_file << M[i];
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, DenseMatrix &M) {
  in_file.setf(ios::skipws);
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/* Define the n x n square tridiagonal matrix class.
   Uses the MV++ MV_Vector_double class. */

/********************************************************
 * Class: TriDiagonalMatrix                             *
 *                                                      *
 * Member functions                                     *
 *      n          -- Return number of rows and columns.*
 *      D          -- Return array of diagonal elements.*
 *      A          -- Return array of elements above    *
 *                    main diagonal.                    *
 *      B          -- Return array of elements below    *
 *                    main diagonal.                    *
 *      allocate   -- Allocate memory for matrix.       *
 *      deallocate -- Deallocate memory for matrix.     *
 *      zero       -- Assign zero matrix.               *
 *      identity   -- Assign identity matrix.           *
 *      trace      -- Return the trace of matrix.       *
 *      transpose  -- Return the transpose of matrix.   *
 *                                                      *
 * Member operators                                     *
 *      M  -- a tridiagonal matrix                      *
 *      RV -- a row vector                              *
 *      CV -- a column vector                           *
 *      a  -- a scalar (double)                         *
 *                                                      *
 * M = M;                                               *
 * RV = M[i];                                           *
 * a = M[i][j];                                         * 
 * M = M + M;                                           *
 * M = M - M;                                           *
 * CV = M * CV;                                         *
 * M = a * M;                                           *
 * M = M * a;                                           *
 * M = M / a;                                           *
 * M = +M;                                              *
 * M = -M;                                              *
 * M += M;                                              *
 * M -= M;                                              *
 * M == M;                                              *
 * M != M;                                              *
 * cout << M; (output function)                         *
 * cin  >> M; (input function)                          *
 *                                                      *
 ********************************************************/
class TriDiagonalMatrix{
  private:
    static RowVector temp_RVec;
  public:
    int                  n;  // size of square tridiagonal matrix.
    MV_Vector_double D,A,B;  // MV++ vectors containing main diagonal
                             // elements and elements above and
                             // below the main diagonal.
                             // Made public so can access them.
			
    /* Creation, copy, and assignment constructors. */
    TriDiagonalMatrix(void) {
       n = 0;
    }

    TriDiagonalMatrix(const TriDiagonalMatrix &M) {
       n = M.n; D = M.D; A = M.A; B = M.B;
    }

    TriDiagonalMatrix(const int nn,
   	              const MV_Vector_double &DD,
	              const MV_Vector_double &AA,
		      const MV_Vector_double &BB) {
       n = nn; D = DD; A = AA; B = BB;
    }

    TriDiagonalMatrix(const int nn,
   	              double *DD,
	              double *AA,
		      double *BB) {
       n = nn; 
       D = MV_Vector_double(DD, nn); 
       A = MV_Vector_double(AA, nn); 
       B = MV_Vector_double(BB, nn);
    }

    /* Destructor. */
    // ~TriDiagonalMatrix(void);
    // Use automatically generated destructor.

    /* Allocate memory for tridiagonal matrix. */
    void allocate(const int N);

    /* Deallocate memory for tridiagonal matrix. */
    void deallocate(void);

    /* Assign the zero matrix. */
    void zero(void);

    /* Assign the identity matrix. */
    void identity(void);

    /* Compute the trace of the matrix. */
    double trace(void);
    double trace(void) const;
    friend double trace(const TriDiagonalMatrix &M);

    /* Assignment operator. */
    // TriDiagonalMatrix operator = (const TriDiagonalMatrix &M);
    // Use automatically generated assignment operator.

    /* Index operator. */
    RowVector &operator[](int index) {
       assert( index >= 0 && index <= n-1 );
       temp_RVec = RowVector(n); temp_RVec.zero();
       if (index > 0 ) temp_RVec(index-1) = B(index-1);
       temp_RVec(index) = D(index);
       if (index < n-1 ) temp_RVec(index+1) = A(index+1);
       return (temp_RVec);
    }
    
    const RowVector &operator[](int index) const {
       assert( index >= 0 && index <= n-1 );
       temp_RVec = RowVector(n); temp_RVec.zero();
       if (index > 0 ) temp_RVec(index-1) = B(index-1);
       temp_RVec(index) = D(index);
       if (index < n-1 ) temp_RVec(index+1) = A(index+1);
       return (temp_RVec);
    }

    /* Unary arithmetic operators. */
    friend TriDiagonalMatrix operator +(const TriDiagonalMatrix &M);
    friend TriDiagonalMatrix operator -(const TriDiagonalMatrix &M);

    /* Binary arithmetic operators. */
    friend TriDiagonalMatrix operator +(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    friend TriDiagonalMatrix operator -(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    friend ColumnVector operator *(const TriDiagonalMatrix &M, const ColumnVector &CVec);
    friend TriDiagonalMatrix operator *(const TriDiagonalMatrix &M, const double &a);
    friend TriDiagonalMatrix operator *(const double &a, const TriDiagonalMatrix &M);
    friend TriDiagonalMatrix operator /(const TriDiagonalMatrix &M, const double &a);

    /* Shortcut arithmetic operators. */
    friend TriDiagonalMatrix &operator +=(TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    friend TriDiagonalMatrix &operator -=(TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    
    /* Relational operators. */
    friend int operator ==(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    friend int operator !=(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const TriDiagonalMatrix &M);
    friend istream &operator >> (istream &in_file, TriDiagonalMatrix &M);

};

/********************************************************
 * TriDiagonalMatrix::allocate -- Allocate memory.      *
 ********************************************************/
inline void TriDiagonalMatrix::allocate(const int N) {
   assert( N >= 1 );
   n = N; D = MV_Vector_double(N, ZERO); 
   A = MV_Vector_double(N, ZERO); B = MV_Vector_double(N, ZERO);
}

/********************************************************
 * TriDiagonalMatrix::deallocate -- Deallocate memory.  *
 ********************************************************/
inline void TriDiagonalMatrix::deallocate(void) {
   n = 0; D = MV_Vector_double(); 
   A = MV_Vector_double(); B = MV_Vector_double();
}

/********************************************************
 * TriDiagonalMatrix::zero -- Assign zero matrix.       *
 ********************************************************/
inline void TriDiagonalMatrix::zero(void) {
   assert( n >= 1 ); int i;
   for ( i = 0; i <= n-1; ++i ) {
      A(i)=ZERO; D(i)=ONE; B(i)=ZERO;
   } /* endfor */
}

/********************************************************
 * TriDiagonalMatrix::indentity -- Set identity matrix. *
 ********************************************************/
inline void TriDiagonalMatrix::identity(void) {
   assert( n >= 1 ); int i;
   for ( i = 0; i <= n-1; ++i ) {
      A(i)=ZERO; D(i)=ONE; B(i)=ZERO;
   } /* endfor */
}

/********************************************************
 * TriDiagonalMatrix::trace -- Trace of a matrix.       *
 ********************************************************/
inline double TriDiagonalMatrix::trace(void) {
   assert( n >= 1 ); int i; double xx; xx= 0;
   for ( i = 0; i <= n-1; ++i ) xx+=D(i);
   return (xx);
}

inline double TriDiagonalMatrix::trace(void) const {
   assert( n >= 1 ); int i; double xx; xx= 0;
   for ( i = 0; i <= n-1; ++i ) xx+=D(i);
   return (xx);
}

inline double trace(const TriDiagonalMatrix &M) {
   assert( M.n >= 1 ); int i; double xx; xx= 0;
   for ( i = 0; i <= M.n-1; ++i ) xx+=M.D(i);
   return (xx);
}

/*!
 * Make matrix M part of the current matrix starting from (StartRow,StartCol) element.
 * This routine is useful is an assembly matrix is formed.
 * \note This routine doesn't check if the current matrix has enough space to 
 *  incorporate matrix M. This property is assumed to exist!
 *
 * \param StartRow the starting matrix row for filling with the elements of M (starts from 0)
 * \param StartCol the starting matrix column for filling with the elements of M (starts from 0)
 */
inline void DenseMatrix::incorporate_matrix(const int & StartRow, const int & StartCol, 
					    const DenseMatrix &M){

  // copy the columns of M into the columns of *this
  // this should run much faster than just accessing each (i,j)
  // element individually 

  int col, row, index;

  for (col = 0; col < M.dim(1); ++col){
    for (row = 0; row < M.lda(); ++row){
      // calculate the index in the current vector
      index = (StartCol + col)*lda() + StartRow + row;
      v_[index] = M.v_[M.lda()*col + row];
    }
  }

}


/********************************************************
 * TriDiagonalMatrix -- Unary arithmetic operators.     *
 ********************************************************/
inline TriDiagonalMatrix operator +(const TriDiagonalMatrix &M) {
   assert( M.n >= 1 );
   return (M);
}

inline TriDiagonalMatrix operator -(const TriDiagonalMatrix &M) {
   assert( M.n >= 1 ); int i; TriDiagonalMatrix M2; M2.allocate(M.n);
   for ( i = 0; i <= M.n-1; ++i ) {
      if (i < M.n-1) M2.A(i) = -M.A(i);  
      M2.D(i) = -M.D(i); 
      if (i > 0) M2.B(i) = -M.B(i);
   } /* endfor */
   return (M2); 
}

/********************************************************
 * TriDiagonalMatrix -- Binary arithmetic operators.    *
 ********************************************************/
inline TriDiagonalMatrix operator +(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   assert( M1.n >= 1 && M1.n == M2.n);
   return (TriDiagonalMatrix(M1.n, M1.D+M2.D, M1.A+M2.A, M1.B+M2.B));
}

inline TriDiagonalMatrix operator -(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   assert( M1.n >= 1 && M1.n == M2.n);
   return (TriDiagonalMatrix(M1.n, M1.D-M2.D, M1.A-M2.A, M1.B-M2.B)); 
}

inline ColumnVector operator *(const TriDiagonalMatrix &M, const ColumnVector &CVec) {
   assert( M.n >= 1 && M.n ==  (int)CVec.size());   
   int i; ColumnVector cv(M.n, ZERO);
   for ( i = 0; i <= M.n-1; ++i ) {
      if (i < M.n-1) cv(i)+=M.A(i)*CVec(i+1);
      cv(i)+=M.D(i)*CVec(i); if (i > 0) cv(i)+=M.B(i)*CVec(i-1);
   } /* endfor */
   return (cv);
}

inline TriDiagonalMatrix operator *(const TriDiagonalMatrix &M, const double &a) {
   assert( M.n >= 1);
   return (TriDiagonalMatrix(M.n, a*M.D, a*M.A, a*M.B));
}

inline TriDiagonalMatrix operator *(const double &a, const TriDiagonalMatrix &M) {
   assert( M.n >= 1);
   return (TriDiagonalMatrix(M.n, a*M.D, a*M.A, a*M.B));
}

inline TriDiagonalMatrix operator /(const TriDiagonalMatrix &M, const double &a) {
   assert( M.n >= 1); double xx; xx = ONE/a;
   return (TriDiagonalMatrix(M.n, xx*M.D, xx*M.A, xx*M.B));
}

/********************************************************
 * TriDiagonalMatrix -- Shortcut arithmetic operators.  *
 ********************************************************/
inline TriDiagonalMatrix &operator +=(TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   assert( M1.n >= 1 && M1.n == M2.n); int i;
   for ( i = 0; i <= M1.n-1; ++i ) {
      if (i < M1.n-1) M1.A(i) += M2.A(i);  
      M1.D(i) += M2.D(i); 
      if (i > 0) M1.B(i) += M2.B(i);
   } /* endfor */
   return (M1);
}

inline TriDiagonalMatrix &operator -=(TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   assert( M1.n >= 1 && M1.n == M2.n); int i;
   for ( i = 0; i <= M1.n-1; ++i ) {
      if (i < M1.n-1) M1.A(i) -= M2.A(i);  
      M1.D(i) -= M2.D(i); 
      if (i > 0) M1.B(i) -= M2.B(i);
   } /* endfor */
   return (M1);
}

/********************************************************
 * TriDiagonalMatrix -- Relational operators.           *
 ********************************************************/
inline int operator ==(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   int i, equal; i = 0; equal = 1;
   if (M1.n != M2.n) equal = 0;
   while (equal) {
      if (i < M1.n-1) { if (M1.A(i) != M2.A(i)) equal = 0; }
      if (M1.D(i) != M2.D(i)) equal = 0;
      if (i > 0) { if (M1.B(i) != M2.B(i)) equal = 0; }
      if (i == M1.n-1) break;
      i = i + 1;
   } /* endwhile */
   return (equal);
}

inline int operator !=(const TriDiagonalMatrix &M1, const TriDiagonalMatrix &M2) {
   int i, not_equal; i = 0; not_equal = 1;
   if (M1.n == M2.n) {
      while (not_equal) {
         if (i < M1.n-1) { if (M1.A(i) == M2.A(i)) not_equal = 0; }
         if (M1.D(i) == M2.D(i)) not_equal = 0;
         if (i > 0) { if (M1.B(i) == M2.B(i)) not_equal = 0; }
         if (i == M1.n-1) break;
         i = i + 1;
      } /* endwhile */
   } /* endif */
   return (not_equal);
}

/********************************************************
 * TriDiagonalMatrix -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const TriDiagonalMatrix &M) {
  int i;
  out_file.setf(ios::scientific);
  for ( i = 0 ; i <= M.n-1; ++i ) {
      out_file << M[i];
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, TriDiagonalMatrix &M) {
  in_file.setf(ios::skipws);
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif /* _MATRIX_INCLUDED  */

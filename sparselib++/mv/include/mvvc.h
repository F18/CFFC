
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//
//      mv_vector_COMPLEX.h       Basic vector class (COMPLEX precision)
//

#ifndef _MV_VECTOR_COMPLEX_H
#define _MV_VECTOR_COMPLEX_H    

// Added by CPTG, Nov. 2003.
#ifndef MVPP_COMPLEX
#ifdef _GNU_GCC_296
#define MVPP_COMPLEX complex<double>
#endif
#ifdef _GNU_GCC_3
#define MVPP_COMPLEX complex<double>
#else
#define MVPP_COMPLEX complex<double>
#endif
#endif

#include <cstdlib>
#include <complex>
#include <iostream>            // for formatted printing of matrices
#ifdef MV_VECTOR_BOUNDS_CHECK  // for vector bounds checking
#   include <cassert>
#endif

using namespace std;

#include "mvvind.h"

// this is really used as a sort of global constant. The reason
// for creating its own type is that so it can be overloaded to perform
// a deep or shallow assignement.  (Any variable of type MV_Vector_::ref_type
// has only one possible value: one.)
//   It is included as a seperate file to avoid multiple definitions.

#include "mvvrf.h"

class MV_Vector_COMPLEX
{                                                                      
    protected:                                                           
           MVPP_COMPLEX *p_;
           unsigned int dim_;
           int ref_;  // 0 or 1; does this own its own memory space?
    public:                                                            


        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    MV_Vector_COMPLEX();                             
    MV_Vector_COMPLEX(unsigned int);                             
    MV_Vector_COMPLEX(unsigned int, const MVPP_COMPLEX&);   
                                                     
                                                    
    MV_Vector_COMPLEX(MVPP_COMPLEX*, unsigned int);      
    MV_Vector_COMPLEX(const MVPP_COMPLEX*, unsigned int);        
    MV_Vector_COMPLEX(const MV_Vector_COMPLEX &); 
    
    // reference of an exisiting data structure
    //
    // note that ref() is initalized with i rather than 1.
    // this is so compilers will not generate a warning that i was
    // not used in the construction.  (MV_Vector::ref_type is an enum that
    // can *only* have the value of 1.
    //
    MV_Vector_COMPLEX(MVPP_COMPLEX* d, unsigned int N, MV_Vector_::ref_type i) :
                            p_(d), dim_(N), ref_(i) {}

    MV_Vector_COMPLEX(const MV_Vector_COMPLEX &V, MV_Vector_::ref_type i)   :
                            p_(V.p_), dim_(V.dim_), ref_(i) {}

    ~MV_Vector_COMPLEX();                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       

    MVPP_COMPLEX&      operator()(unsigned int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  MVPP_COMPLEX&       operator()(unsigned int i) const 
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }

    MVPP_COMPLEX&      operator[](unsigned int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  MVPP_COMPLEX&       operator[](unsigned int i) const 
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }


    MV_Vector_COMPLEX operator()(const MV_VecIndex &I) ;
    MV_Vector_COMPLEX operator()(void);
    const MV_Vector_COMPLEX operator()(void) const;
    const MV_Vector_COMPLEX operator()(const MV_VecIndex &I) const;

    inline unsigned int             size() const { return dim_;}
    inline unsigned int             dim() const { return dim_;}
    inline int                      ref() const { return  ref_;}
    inline int                      null() const {return dim_== 0;}
            //
            // Create a new *uninitalized* vector of size N
            MV_Vector_COMPLEX & newsize(unsigned int );
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    MV_Vector_COMPLEX & operator=(const MV_Vector_COMPLEX&);
    MV_Vector_COMPLEX & operator=(const MVPP_COMPLEX&);


    friend ostream& operator<<(ostream &s, const MV_Vector_COMPLEX &A);

};                                                                     



#endif  

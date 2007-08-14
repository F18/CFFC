/* Euler3DState.h:  Header file defining 3D Euler Solution State Classes. */

#ifndef _EULER3D_STATE_INCLUDED
#define _EULER3D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, 3D vector, and gas constant header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/* Define the classes. */

#define	NUM_VAR_EULER3D    5

class Euler3D_cState;

/********************************************************
 * Class: Euler3D_pState                                *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1                           *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      E       -- Return total energy.                 *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      dv      -- Return momentum.                     *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          *
 *      U       -- Return conserved solution state.     *
 *      Fx      -- Return x-direction solution flux.    *
 *      Fy      -- Return y-direction solution flux.    *
 *      Fz      -- Return z-direction solution flux.    *
 *      Fn      -- Return n-direction solution flux.    *
 *      lambda  -- Return eigenvalue.                   *
 *      rp      -- Return primitive right eigenvector.  *
 *      rc      -- Return conserved right eigenvector.  *
 *      lp      -- Return primitive left eigenvector.   *
 *                                                      *
 * Member operators                                     *
 *      W -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * W = W;                                               *
 * c = W[i];                                            *
 * W = W + W;                                           *
 * W = W - W;                                           *
 * c = W * W; (inner product)                           *
 * W = c * W;                                           *
 * W = W * c;                                           *
 * W = W / c;                                           *
 * W = +W;                                              *
 * W = -W;                                              *
 * W += W;                                              *
 * W -= W;                                              *
 * W == W;                                              *
 * W != W;                                              *
 * cout << W; (output function)                         *
 * cin  >> W; (input function)                          *
 *                                                      *
 ********************************************************/
class Euler3D_pState{
  private:

   enum alloc_type_t {
      not_allocated,
      alloc_new,
      alloc_raw
   };
   alloc_type_t allocType;

  public:
   
   double       *solnvec;    //Solution vector
   static double   g;   // Specific heat ratio.
   static double gm1;   // g-1
   static double gm1i;  // 1/(g-1)
   static double   R;   // Gas constant.
   // Made public so can access them.
   /* Initialization of the solnvec */
   /* make the pointer (solnvec) pointing to the first element */
   /* of the solution block memory for primitive solution vector*/
   void Initialization(double *first_element_in_W_memory ){ 
      solnvec = first_element_in_W_memory;
      allocType = alloc_raw;
   } 
   /* creation constructor */

   // Default constructor
   Euler3D_pState(void) : allocType(not_allocated), solnvec(0) { }

   // Construct using new or address
   Euler3D_pState(double *const address) {
      if ( address ) {
         // solnvec is allocated on given address
         solnvec = address;
         allocType = alloc_raw;
      }
      else {
         //solnvec is allocated using new 
         solnvec = new double[NUM_VAR_EULER3D];
         allocType = alloc_new;
         //assignment
         d() = DENSITY_STDATM;
         vx() = ZERO;  vy() = ZERO;
         vz() = ZERO;
         p() = PRESSURE_STDATM;
      }
   }

   // Copy constructor (uses new)
   Euler3D_pState(const Euler3D_pState &W) {
      solnvec = new double[NUM_VAR_EULER3D];
      allocType = alloc_new;
      for ( int i = 0 ; i < NUM_VAR_EULER3D ; ++i )
         solnvec[i] = W.solnvec[i];
   }

   // Assignment operator
   Euler3D_pState& operator=(const Euler3D_pState &W) {
      if ( this !=&W ) {
         assert(W.allocType != not_allocated);
         if ( allocType == not_allocated ) {
            solnvec = new double[NUM_VAR_EULER3D];
            allocType = alloc_new;
         }
       
         for ( int i = 0 ; i != NUM_VAR_EULER3D ; ++i )
            solnvec[i] = W.solnvec[i];
        
      }

      return (*this);
   }


   // Destructor
   ~Euler3D_pState() {
      if ( allocType ==  alloc_new ) {
         delete[] solnvec;
         allocType = not_allocated;
      }
   }
   
    
   Euler3D_pState(const Euler3D_cState &U);
   
   Euler3D_pState(const double &rho,
                  const Vector3D &V,
                  const double &pre) {
      
      solnvec = new double[NUM_VAR_EULER3D];
      allocType = alloc_new;
      d() = rho;  vx() = V.x;
      vy() = V.y; vz() = V.z;
      p() = pre;
   }
   
   Euler3D_pState(const double &rho,
                  const double &Vx,
                  const double &Vy,
                  const double &Vz,
                  const double &pre) {
      
      solnvec = new double[NUM_VAR_EULER3D];
      allocType = alloc_new;
      
       d()  = rho;  vx() = Vx;
       vy() = Vy;   vz() = Vz; 
       p() = pre;
      
   }
   
   /* Vacuum operator. */
   void Vacuum() {
      for( int i= 0; i!= NUM_VAR_EULER3D; ++i)
         solnvec[i] =  ZERO; 
      
   }
   /* Destructor. */
   // ~Euler3D_pState(void);
   // Use automatically generated destructor.
   
   /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

   /* Total velocity */
    double uo(void) const;
    
    /* Temperature. */
    double T(void);
    double T(void) const;
    
    /* Specific internal energy. */
    double e(void);
    double e(void) const;

    /* Total energy. */
    double E(void);
    double E(void) const;
    
    /* Specific enthalpy. */
    double h(void);
    double h(void) const;

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Mach number. */
    double M(void);
    double M(void) const;
    
    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Momentum. */
    Vector3D dv(void);
    Vector3D dv(void) const;
    double dv(const Vector3D &n);
    double dv(const Vector3D &n) const;

    /* Stagnation temperature. */
    double To(void);
    double To(void) const;

    /* Stagnation pressure. */
    double po(void);
    double po(void) const;

    /* Stagnation sound speed. */
    double ao(void);
    double ao(void) const;

    /* Stagnation enthalpy. */
    double ho(void);
    double ho(void) const;

    double &d(void){ return solnvec[0]; }
    double d(void) const{ return solnvec[0]; }
    double &vx(void){ return solnvec[1]; }
    double vx(void) const{ return solnvec[1]; }
    double &vy(void){ return solnvec[2]; }
    double vy(void) const{ return solnvec[2]; }
    double &vz(void){ return solnvec[3]; }
    double vz(void) const{ return solnvec[3]; }
    double &p(void){ return solnvec[4]; }
    double p(void) const{ return solnvec[4]; }
  
    /* Conserved solution state. */
    Euler3D_cState U(void);
    Euler3D_cState U(void) const;
    Euler3D_cState U(const Euler3D_pState &W);
    friend Euler3D_cState U(const Euler3D_pState &W);
    
    /* Assignment operator. */
    // Euler3D_pState operator = (const Euler3D_pState &W);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
       assert( index >= 1 && index <= NUM_VAR_EULER3D );
       
       return solnvec[index-1];
    }
    
    double operator[](int index) const {
       assert( index >= 1 && index <= NUM_VAR_EULER3D );
       return solnvec[index-1];
    }

    /* Eigenvalue(s) (x-direction). */
    Euler3D_pState lambda(void);
    Euler3D_pState lambda(void) const;
    Euler3D_pState lambda(const Euler3D_pState &W);
    friend Euler3D_pState lambda(const Euler3D_pState &W);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const Euler3D_pState &W, int index);

    Euler3D_pState lambda_x(void);
    Euler3D_pState lambda_x(void) const;
    Euler3D_pState lambda_x(const Euler3D_pState &W);
    friend Euler3D_pState lambda_x(const Euler3D_pState &W);
    double lambda_x(int index);
    double lambda_x(int index) const;
    friend double lambda_x(const Euler3D_pState &W, int index);

    /* Binary arithmetic operators. */
    friend Euler3D_pState operator +(const Euler3D_pState &W1, 
                                     const Euler3D_pState &W2);
    friend Euler3D_pState operator -(const Euler3D_pState &W1, 
                                     const Euler3D_pState &W2);
    friend double operator *(const Euler3D_pState &W1, 
                             const Euler3D_pState &W2);
    friend Euler3D_pState operator *(const Euler3D_pState &W, const double &a);
    friend Euler3D_pState operator *(const double &a, const Euler3D_pState &W);
    friend Euler3D_pState operator /(const Euler3D_pState &W, const double &a);

    /* Unary arithmetic operators. */
    friend Euler3D_pState operator +(const Euler3D_pState &W);
    friend Euler3D_pState operator -(const Euler3D_pState &W);

    /* Shortcut arithmetic operators. */
    friend Euler3D_pState &operator +=(Euler3D_pState &W1, 
                                       const Euler3D_pState &W2);
    friend Euler3D_pState &operator -=(Euler3D_pState &W1, 
                                       const Euler3D_pState &W2);
    
    /* Relational operators. */
    friend int operator ==(const Euler3D_pState &W1, const Euler3D_pState &W2);
    friend int operator !=(const Euler3D_pState &W1, const Euler3D_pState &W2);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Euler3D_pState &W);
    friend istream &operator >> (istream &in_file,  Euler3D_pState &W);
    friend Euler3D_pState operator ^(const Euler3D_pState &W1, 
                                     const Euler3D_pState &W2);
    


    Euler3D_cState F(void);
    Euler3D_cState F(void) const;
    Euler3D_cState F(const Euler3D_pState &W);
    friend Euler3D_cState F(const Euler3D_pState &W);

    Euler3D_cState Fx(void);
    Euler3D_cState Fx(void) const;
    Euler3D_cState Fx(const Euler3D_pState &W);
    friend Euler3D_cState Fx(const Euler3D_pState &W);
    
};

/********************************************************
 * Class: Euler3D_cState                                *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      E       -- Return total energy.                 *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1.                          *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          *
 *      W       -- Return primitive solution state.     *
 *                                                      *
 * Member operators                                     *
 *      U -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * U = U;                                               *
 * c = U[i];                                            *
 * U = U + U;                                           *
 * U = U - U;                                           *
 * c = U * U; (inner product)                           *
 * U = c * U;                                           *
 * U = U * c;                                           *
 * U = U / c;                                           *
 * U = +U;                                              *
 * U = -U;                                              *
 * U += U;                                              *
 * U -= U;                                              *
 * U == U;                                              *
 * U != U;                                              *
 * cout << U; (output function)                         *
 * cin  >> U; (input function)                          *
 *                                                      *
 ********************************************************/
class Euler3D_cState{
  private:
  public:
    double          d;   // Density.
    Vector3D       dv;   // Momentum.
    double          E;   // Total Energy.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
    static double   R;   // Gas constant.
	                 // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Euler3D_cState(void) {
       d = DENSITY_STDATM; dv.zero(); E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
    }

    Euler3D_cState(const Euler3D_cState &U) {
       d = U.d; dv = U.dv; E = U.E;
    }

    Euler3D_cState(const Euler3D_pState &W);

    Euler3D_cState(const double &rho,
	           const Vector3D &rhoV,
	           const double &Etotal) {
       d = rho; dv = rhoV; E = Etotal;
    }

    Euler3D_cState(const double &rho,
	           const double &rhovx,
	           const double &rhovy,
	           const double &rhovz,
	           const double &Etotal) {
       d = rho; dv.x = rhovx; dv.y = rhovy; dv.z = rhovz; E = Etotal;
    }
    
    /* Destructor. */
    // ~Euler3D_cState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Flow velocity. */
    Vector3D v(void);
    Vector3D v(void) const;
    double v(const Vector3D &n);
    double v(const Vector3D &n) const;
    
    /* Pressure. */
    double p(void);
    double p(void) const;

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Specific internal energy. */
    double e(void);
    double e(void) const;

    /* Specific enthalpy. */
    double h(void);
    double h(void) const;

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Mach number. */
    double M(void);
    double M(void) const;

    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Stagnation temperature. */
    double To(void);
    double To(void) const;

    /* Stagnation pressure. */
    double po(void);
    double po(void) const;

    /* Stagnation sound speed. */
    double ao(void);
    double ao(void) const;

    /* Stagnation enthalpy. */
    double ho(void);
    double ho(void) const;

    

    /* Primitive solution state. */
    Euler3D_pState W(void);
    Euler3D_pState W(void) const;
    Euler3D_pState W(const Euler3D_cState &U);
    friend Euler3D_pState W(const Euler3D_cState &U);
    
    /* Assignment operator. */
    // Euler3D_cState operator = (const Euler3D_cState &U);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_EULER3D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (dv.z);
        case 5 :
	  return (E);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_EULER3D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (dv.z);
        case 5 :
	  return (E);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend Euler3D_cState operator +(const Euler3D_cState &U1, const Euler3D_cState &U2);
    friend Euler3D_cState operator -(const Euler3D_cState &U1, const Euler3D_cState &U2);
    friend double operator *(const Euler3D_cState &U1, const Euler3D_cState &U2);
    friend Euler3D_cState operator *(const Euler3D_cState &U, const double &a);
    friend Euler3D_cState operator *(const double &a, const Euler3D_cState &U);
    friend Euler3D_cState operator /(const Euler3D_cState &U, const double &a);

    /* Unary arithmetic operators. */
    friend Euler3D_cState operator +(const Euler3D_cState &U);
    friend Euler3D_cState operator -(const Euler3D_cState &U);

    /* Shortcut arithmetic operators. */
    friend Euler3D_cState &operator +=(Euler3D_cState &U1, const Euler3D_cState &U2);
    friend Euler3D_cState &operator -=(Euler3D_cState &U1, const Euler3D_cState &U2);
    
    /* Relational operators. */
    friend int operator ==(const Euler3D_cState &U1, const Euler3D_cState &U2);
    friend int operator !=(const Euler3D_cState &U1, const Euler3D_cState &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Euler3D_cState &U);
    friend istream &operator >> (istream &in_file,  Euler3D_cState &U);
    friend Euler3D_cState operator ^(const Euler3D_cState &U1, const Euler3D_cState &U2);
    
    Euler3D_cState Fx(void);
    Euler3D_cState Fx(void) const;
    Euler3D_cState Fx(const Euler3D_cState &U);
    friend Euler3D_cState Fx(const Euler3D_cState &U);
       
};

/********************************************************
 * Euler3D_pState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler3D_pState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler3D_pState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } else if (strcmp(string_ptr, "A") == 0) {
     g = GAMMA_A;
     R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
   } else if (strcmp(string_ptr, "CO") == 0) {
     g = GAMMA_CO;
     R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
   } else if (strcmp(string_ptr, "CO2") == 0) {
     g = GAMMA_CO2;
     R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
   } else if (strcmp(string_ptr, "CH4") == 0) {
     g = GAMMA_CH4;
     R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
   } else if (strcmp(string_ptr, "H") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
   } else if (strcmp(string_ptr, "H2") == 0) {
     g = GAMMA_H2;
     R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
   } else if (strcmp(string_ptr, "HE") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
   } else if (strcmp(string_ptr, "H2O") == 0) {
     g = GAMMA_H2O;
     R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
   } else if (strcmp(string_ptr, "N2") == 0) {
     g = GAMMA_N2;
     R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
   } else if (strcmp(string_ptr, "O") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
   } else if (strcmp(string_ptr, "O2") == 0) {
     g = GAMMA_O2;
     R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
   } else if (strcmp(string_ptr, "e") == 0) {
     g = GAMMA_e;
     R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler3D_pState::uo -- Total velocity                 *
 ********************************************************/
   /* Total velocity */
inline double Euler3D_pState::uo(void) const{
   return (sqrt(vx()*vx() +vy()*vy() +vz()*vz()));   
   
}

/********************************************************
 * Euler3D_pState::T -- Temperature.                    *
 ********************************************************/
inline double Euler3D_pState::T(void) {
    return (p()/(d()*R));
}

inline double Euler3D_pState::T(void) const {
    return (p()/(d()*R));
}

/********************************************************
 * Euler3D_pState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler3D_pState::e(void) {
    return (p()/(gm1*d()));
}

inline double Euler3D_pState::e(void) const {
    return (p()/(gm1*d()));
}

/********************************************************
 * Euler3D_pState::E -- Total energy.                   *
 ********************************************************/
inline double Euler3D_pState::E(void) {
   return (p()*gm1i + HALF*d()*sqr(uo()));
}

inline double Euler3D_pState::E(void) const {
   return (p()*gm1i + HALF*d()*sqr(uo()));
}

/********************************************************
 * Euler3D_pState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler3D_pState::h(void) {
    return (g*p()/(gm1*d()) + HALF*sqr(uo()));
}

inline double Euler3D_pState::h(void) const {
    return (g*p()/(gm1*d()) + HALF*(sqr(uo())));
}

/********************************************************
 * Euler3D_pState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler3D_pState::H(void) {
    return (g*gm1i*p() + HALF*d()*sqr(uo()));
}

inline double Euler3D_pState::H(void) const {
    return (g*gm1i*p() + HALF*d()*sqr(uo()));
}

/********************************************************
 * Euler3D_pState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler3D_pState::a(void) {
    return (sqrt(g*p()/d()));
}

inline double Euler3D_pState::a(void) const {
    return (sqrt(g*p()/d()));
}

/********************************************************
 * Euler3D_pState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler3D_pState::a2(void) {
    return (g*p()/d());
}

inline double Euler3D_pState::a2(void) const {
    return (g*p()/d());
}

/********************************************************
 * Euler3D_pState::M -- Mach number.                    *
 ********************************************************/
inline double Euler3D_pState::M(void) {
    return (uo()/sqrt(g*p()/d()));
}

inline double Euler3D_pState::M(void) const {
    return (uo()/sqrt(g*p()/d()));
}

/********************************************************
 * Euler3D_pState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler3D_pState::s(void) {
    return (R*gm1i*log(p()/pow(d(), g)));
}

inline double Euler3D_pState::s(void) const {
    return (R*gm1i*log(p()/pow(d(), g)));
}

/********************************************************
 * Euler3D_pState::dv -- Momentum.                      *
 ********************************************************/
inline Vector3D Euler3D_pState::dv(void) {
   Vector3D V;
   V.x = vx(); V.y = vy(); V.z = vz();
   
   return (d()*V);
}

inline Vector3D Euler3D_pState::dv(void) const {
   Vector3D V;
   V.x = vx(); V.y = vy(); V.z = vz();
   return (d()*V);
   
}

inline double Euler3D_pState::dv(const Vector3D &n) {
   Vector3D V;
   V.x = vx(); V.y = vy(); V.z = vz();
   return (d()*(V*n));
   
}

inline double Euler3D_pState::dv(const Vector3D &n) const {
   Vector3D V;
   V.x = vx(); V.y = vy(); V.z = vz();
   return (d()*(V*n));
   
}

/********************************************************
 * Euler3D_pState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler3D_pState::To(void) {
   return ((p()/(d()*R))*(ONE+HALF*gm1*sqr(uo())/(g*p()/d())));
}

inline double Euler3D_pState::To(void) const {
   return ((p()/(d()*R))*(ONE+HALF*gm1*sqr(uo())/(g*p()/d())));
}

/********************************************************
 * Euler3D_pState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler3D_pState::po(void) {
   return (p()*pow(ONE+HALF*gm1*sqr(uo())/(g*p()/d()), g*gm1i));
}

inline double Euler3D_pState::po(void) const {
   return (p()*pow(ONE+HALF*gm1*sqr(uo())/(g*p()/d()), g*gm1i));
}

/********************************************************
 * Euler3D_pState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler3D_pState::ao(void) {
   return (sqrt((g*p()/d())*(ONE+HALF*gm1*sqr(uo())/(g*p()/d()))));
}

inline double Euler3D_pState::ao(void) const {
   return (sqrt((g*p()/d())*(ONE+HALF*gm1*sqr(uo())/(g*p()/d()))));
}

/********************************************************
 * Euler3D_pState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler3D_pState::ho(void) {
   return ((g*p()/(gm1*d()) + HALF*sqr(uo()))*
           (ONE+HALF*gm1*sqr(uo())/(g*p()/d())));
}

inline double Euler3D_pState::ho(void) const {
   return ((g*p()/(gm1*d()) + HALF*sqr(uo()))
           *(ONE+HALF*gm1*sqr(uo())/(g*p()/d())));
}

/********************************************************
 * Euler3D_pState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler3D_pState operator +(const Euler3D_pState &W1, 
                                 const Euler3D_pState &W2) {
   return (Euler3D_pState(W1.d()+W2.d(), W1.vx()+W2.vx(),
                          W1.vy()+W2.vy(), W1.vz()+W2.vz(),
                          W1.p()+W2.p()));
}

inline Euler3D_pState operator -(const Euler3D_pState &W1, 
                                 const Euler3D_pState &W2) {
   return (Euler3D_pState(W1.d()-W2.d(),W1.vx()-W2.vx(),
                          W1.vy()-W2.vy(), W1.vz()-W2.vz(),
                          W1.p()-W2.p()));
}

// Inner product operator.
inline double operator *(const Euler3D_pState &W1, const Euler3D_pState &W2) {
   return (W1.d()*W2.d()+W1.vx()*W2.vx()+
           W1.vy()*W2.vy() + W1.vz()*W2.vz(),
           +W1.p()*W2.p());
}

inline Euler3D_pState operator *(const Euler3D_pState &W, const double &a) {
   return (Euler3D_pState(a*W.d(),a*W.vx(),a*W.vy(), a*W.vz(), a*W.p()));
}

inline Euler3D_pState operator *(const double &a, const Euler3D_pState &W) {
  return (Euler3D_pState(a*W.d(),a*W.vx(),a*W.vy(), a*W.vz(),a*W.p()));
}

inline Euler3D_pState operator /(const Euler3D_pState &W, const double &a) {
  return (Euler3D_pState(W.d()/a, W.vx()/a, W.vy()/a, W.vz()/a,W.p()/a));
}
inline Euler3D_pState operator ^(const Euler3D_pState &W1, 
                                 const Euler3D_pState &W2) {
   return (Euler3D_pState(W1.d()*W2.d(),W1.vx()*W2.vx(),
                          W1.vy()*W2.vy(), 
                          W1.vz()*W2.vz(), 
                          W1.p()*W2.p()));
}
/********************************************************
 * Euler3D_pState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler3D_pState operator +(const Euler3D_pState &W) {
   return (Euler3D_pState(W.d(),W.vx(),W.vy(), W.vz(),W.p()));
}

inline Euler3D_pState operator -(const Euler3D_pState &W) {
   return (Euler3D_pState(-W.d(),-W.vx(), -W.vy(), -W.vz(), -W.p()));
}

/********************************************************
 * Euler3D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler3D_pState &operator +=(Euler3D_pState &W1,
                                   const Euler3D_pState &W2) {

   for ( int i = 0 ; i != NUM_VAR_EULER3D ; ++i ) 
      W1.solnvec[i] += W2.solnvec[i];
   
   return (W1);
}

inline Euler3D_pState &operator -=(Euler3D_pState &W1, 
                                   const Euler3D_pState &W2) {
   
   for ( int i = 0 ; i != NUM_VAR_EULER3D ; ++i ) 
      W1.solnvec[i] -= W2.solnvec[i];
   
   return (W1);
}

/************************************************************
 * Euler3D_pState::lambda -- Eigenvalue(s) (x-direction).   *
 ************************************************************/
inline Euler3D_pState Euler3D_pState::lambda(void) {
   double c = a();
   return (Euler3D_pState(vx() - c, vx(), vx(), vx(), vx() + c));
}

inline Euler3D_pState Euler3D_pState::lambda(void) const {
   double c = a();
   return (Euler3D_pState(vx() - c, vx(), vx(), vx(), vx() + c));
}

inline Euler3D_pState Euler3D_pState::lambda(const Euler3D_pState &W) {
   double c = W.a();
   return (Euler3D_pState(W.vx() - c, W.vx(), W.vx(), W.vx(), W.vx() + c));
}

inline Euler3D_pState lambda(const Euler3D_pState &W) {
   double c = W.a();
   return (Euler3D_pState(W.vx() - c, W.vx(), W.vx(), W.vx(), W.vx() + c));
}

inline double Euler3D_pState::lambda(int index) {
   assert( index >= 1 && index <= NUM_VAR_EULER3D );
   switch(index) {
   case 1 :
      return (vx()-a());
   case 2 :
      return (vx());
   case 3 :
      return (vx());
   case 4 :
      return (vx());
   case 5 :
      return (vx()+a());
   default:
      return (vx());
   };
}

inline double Euler3D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_EULER3D );
  switch(index) {
  case 1 :
     return (vx()-a());
  case 2 :
     return (vx());
  case 3 :
     return (vx());
  case 4 :
     return (vx());
  case 5 :
     return (vx()+a());
  default:
     return (vx());
  };
}

inline double lambda(const Euler3D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER3D );
  switch(index) {
  case 1 :
     return (W.vx()-W.a());
  case 2 :
     return (W.vx());
  case 3 :
     return (W.vx());
  case 4 :
     return (W.vx());
  case 5 :
     return (W.vx()+W.a());
  default:
     return (W.vx());
  };
}

/************************************************************
 * Euler3D_pState::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline Euler3D_pState Euler3D_pState::lambda_x(void) {
  double c = a();
  return (Euler3D_pState(vx() - c, vx(), vx(), vx(), vx() + c));
}

inline Euler3D_pState Euler3D_pState::lambda_x(void) const {
  double c = a();
  return (Euler3D_pState(vx() - c, vx(), vx(), vx(), vx() + c));
}

inline Euler3D_pState Euler3D_pState::lambda_x(const Euler3D_pState &W) {
  double c = W.a();
  return (Euler3D_pState(W.vx() - c, W.vx(), W.vx(), W.vx(), W.vx() + c));
}

inline Euler3D_pState lambda_x(const Euler3D_pState &W) {
  double c = W.a();
  return (Euler3D_pState(W.vx() - c, W.vx(), W.vx(), W.vx(), W.vx() + c));
}

inline double Euler3D_pState::lambda_x(int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER3D );
  switch(index) {
  case 1 :
     return (vx()-a());
  case 2 :
     return (vx());
  case 3 :
     return (vx());
  case 4 :
     return (vx());
  case 5 :
     return (vx()+a());
  default:
     return (vx());
  };
}

inline double Euler3D_pState::lambda_x(int index) const {
   assert( index >= 1 && index <= NUM_VAR_EULER3D );
   switch(index) {
   case 1 :
      return (vx()-a());
   case 2 :
      return (vx());
   case 3 :
      return (vx());
   case 4 :
      return (vx());
   case 5 :
      return (vx()+a());
   default:
      return (vx());
   };
}

inline double lambda_x(const Euler3D_pState &W, int index) {
   assert( index >= 1 && index <= NUM_VAR_EULER3D);
   switch(index) {
   case 1 :
      return (W.vx()-W.a());
   case 2 :
      return (W.vx());
   case 3 :
      return (W.vx());
   case 4 :
      return (W.vx());
   case 5 :
      return (W.vx()+W.a());
   default:
      return (W.vx());
   };
}

/********************************************************
 * Euler3D_pState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler3D_pState &W1, const Euler3D_pState &W2) {
   return (W1.d() == W2.d() && W1.vx() == W2.vx() && W1.vy() == W2.vy()
           && W1.vz() == W2.vz() && W1.p() == W2.p());
}

inline int operator !=(const Euler3D_pState &W1, const Euler3D_pState &W2) {
   return (W1.d() != W2.d() || W1.vx() != W2.vx() || W1.vy() != W2.vy()
           || W1.vz() != W2.vz() || W1.p() != W2.p());
}

/********************************************************
 * Euler3D_pState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler3D_pState &W) {
   out_file.setf(ios::scientific);
   out_file << " " << W.d()<< " " << W.vx() << " " << W.vy() 
            << " " << W.vz() << " " << W.p();
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, Euler3D_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.d() >> W.vx() >> W.vy() >> W.vz() 
           >> W.p();
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/********************************************************
 * Euler3D_cState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler3D_cState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler3D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } else if (strcmp(string_ptr, "A") == 0) {
     g = GAMMA_A;
     R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
   } else if (strcmp(string_ptr, "CO") == 0) {
     g = GAMMA_CO;
     R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
   } else if (strcmp(string_ptr, "CO2") == 0) {
     g = GAMMA_CO2;
     R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
   } else if (strcmp(string_ptr, "CH4") == 0) {
     g = GAMMA_CH4;
     R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
   } else if (strcmp(string_ptr, "H") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
   } else if (strcmp(string_ptr, "H2") == 0) {
     g = GAMMA_H2;
     R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
   } else if (strcmp(string_ptr, "HE") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
   } else if (strcmp(string_ptr, "H2O") == 0) {
     g = GAMMA_H2O;
     R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
   } else if (strcmp(string_ptr, "N2") == 0) {
     g = GAMMA_N2;
     R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
   } else if (strcmp(string_ptr, "O") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
   } else if (strcmp(string_ptr, "O2") == 0) {
     g = GAMMA_O2;
     R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
   } else if (strcmp(string_ptr, "e") == 0) {
     g = GAMMA_e;
     R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler3D_cState::v -- Flow velocity.                  *
 ********************************************************/
inline Vector3D Euler3D_cState::v(void) {
    return (dv/d);
}

inline Vector3D Euler3D_cState::v(void) const {
    return (dv/d);
}

inline double Euler3D_cState::v(const Vector3D &n) {
    return ((dv*n)/d);
}

inline double Euler3D_cState::v(const Vector3D &n) const {
    return ((dv*n)/d);
}

/********************************************************
 * Euler3D_cState::p -- Pressure.                       *
 ********************************************************/
inline double Euler3D_cState::p(void) {
    return (gm1*(E - HALF*dv.sqr()/d));
}

inline double Euler3D_cState::p(void) const {
    return (gm1*(E - HALF*dv.sqr()/d));
}

/********************************************************
 * Euler3D_cState::T -- Temperature.                    *
 ********************************************************/
inline double Euler3D_cState::T(void) {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

inline double Euler3D_cState::T(void) const {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

/********************************************************
 * Euler3D_cState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler3D_cState::e(void) {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

inline double Euler3D_cState::e(void) const {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

/********************************************************
 * Euler3D_cState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler3D_cState::h(void) {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

inline double Euler3D_cState::h(void) const {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

/********************************************************
 * Euler3D_cState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler3D_cState::H(void) {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

inline double Euler3D_cState::H(void) const {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

/********************************************************
 * Euler3D_cState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler3D_cState::a(void) {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

inline double Euler3D_cState::a(void) const {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

/********************************************************
 * Euler3D_cState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler3D_cState::a2(void) {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

inline double Euler3D_cState::a2(void) const {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

/********************************************************
 * Euler3D_cState::M -- Mach number.                    *
 ********************************************************/
inline double Euler3D_cState::M(void) {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler3D_cState::M(void) const {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler3D_cState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler3D_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

inline double Euler3D_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

/********************************************************
 * Euler3D_cState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler3D_cState::To(void) {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler3D_cState::To(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler3D_cState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler3D_cState::po(void) {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*
                                       (E/d -HALF*dv.sqr()/sqr(d))), g*gm1i));
}

inline double Euler3D_cState::po(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*
                                       (E/d - HALF*dv.sqr()/sqr(d))), g*gm1i));
}

/********************************************************
 * Euler3D_cState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler3D_cState::ao(void) {
   return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
                (ONE+HALF*gm1*dv.sqr()/
                 (d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

inline double Euler3D_cState::ao(void) const {
   return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
                (ONE+HALF*gm1*dv.sqr()/
                 (d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

/********************************************************
 * Euler3D_cState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler3D_cState::ho(void) {
   return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
           (ONE+HALF*gm1*dv.sqr()/
            (d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler3D_cState::ho(void) const {
    return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
	    (ONE+HALF*gm1*dv.sqr()/
             (d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler3D_cState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler3D_cState operator +(const Euler3D_cState &U1, 
                                 const Euler3D_cState &U2) {
   return (Euler3D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E));
}

inline Euler3D_cState operator -(const Euler3D_cState &U1, 
                                 const Euler3D_cState &U2) {
   return (Euler3D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E));
}

// Inner product operator.
inline double operator *(const Euler3D_cState &U1, 
                         const Euler3D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E*U2.E);
}

inline Euler3D_cState operator *(const Euler3D_cState &U, const double &a) {
   return (Euler3D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler3D_cState operator *(const double &a, const Euler3D_cState &U) {
   return (Euler3D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler3D_cState operator /(const Euler3D_cState &U, const double &a) {
   return (Euler3D_cState(U.d/a,U.dv/a,U.E/a));
}

inline Euler3D_cState operator ^(const Euler3D_cState &U1, 
                                 const Euler3D_cState &U2) {
   return (Euler3D_cState(U1.d*U2.d,U1.dv.x*U2.dv.x,
                          U1.dv.y*U2.dv.y,U1.dv.z*U2.dv.z, 
                          U1.E*U2.E));
}

/********************************************************
 * Euler3D_cState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler3D_cState operator +(const Euler3D_cState &U) {
   return (Euler3D_cState(U.d,U.dv,U.E));
}

inline Euler3D_cState operator -(const Euler3D_cState &U) {
   return (Euler3D_cState(-U.d,-U.dv,-U.E));
}

/********************************************************
 * Euler3D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler3D_cState &operator +=(Euler3D_cState &U1, const Euler3D_cState &U2) {
   U1.d += U2.d;
   U1.dv += U2.dv;
   U1.E += U2.E;
   return (U1);
}

inline Euler3D_cState &operator -=(Euler3D_cState &U1, const Euler3D_cState &U2) {
   U1.d -= U2.d;
   U1.dv -= U2.dv;
   U1.E -= U2.E;
   return (U1);
}

/********************************************************
 * Euler3D_cState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler3D_cState &U1, const Euler3D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E);
}

inline int operator !=(const Euler3D_cState &U1, const Euler3D_cState &U2) {
   return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E);
}

/********************************************************
 * Euler3D_cState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler3D_cState &U) {
  
   out_file.setf(ios::scientific);
   out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " 
            << U.dv.z << " " << U.E;
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, Euler3D_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.d >> U.dv.x >> U.dv.y >> U.dv.z >> U.E;
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/********************************************************
 * Euler3D_pState::Euler3D_pState -- Constructor.       *
 ********************************************************/
inline Euler3D_pState::Euler3D_pState(const Euler3D_cState &U) {
   solnvec = new double[NUM_VAR_EULER3D];
   allocType = alloc_new;
   d() = U.d; vx() = U.v().x;
   vy()= U.v().y; vz() = U.v().z;
   p() = U.p();
}

/********************************************************
 * Euler3D_pState::U -- Conserved solution state.       *
 ********************************************************/
inline Euler3D_cState Euler3D_pState::U(void) {
   return (Euler3D_cState(d(), dv(), E()));
}

inline Euler3D_cState Euler3D_pState::U(void) const {
   return (Euler3D_cState(d(), dv(), E()));
}

inline Euler3D_cState Euler3D_pState::U(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d(), W.d()*W.vx(), W.d()*W.vy(),
			  W.d()*W.vz(), W.E()));
}

inline Euler3D_cState U(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d(), W.dv(), W.E()));
}

/********************************************************
 * Euler3D_cState::Euler3D_cState -- Constructor.       *
 ********************************************************/
inline Euler3D_cState::Euler3D_cState(const Euler3D_pState &W) {
   d = W.d(); dv = W.dv(); E = W.E();
}

/********************************************************
 * Euler3D_cState::W -- Primitive solution state.       *
 ********************************************************/
inline Euler3D_pState Euler3D_cState::W(void) {
   return (Euler3D_pState(d, v(), p()));
}

inline Euler3D_pState Euler3D_cState::W(void) const {
  return (Euler3D_pState(d, v(), p()));
}

inline Euler3D_pState Euler3D_cState::W(const Euler3D_cState &U) {
  return (Euler3D_pState(U.d, U.v(), U.p()));
}

inline Euler3D_pState W(const Euler3D_cState &U) {
   return (Euler3D_pState(U.d, U.v(), U.p()));
}

/********************************************************
 * Euler3D_pState::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline Euler3D_cState Euler3D_pState::F(void) {
   return (Euler3D_cState(d()*vx(), d()*sqr(vx()) + p(),
                          d()*vx()*vy(), d()*vx()*vz(), 
                          vx()*H()));
}

inline Euler3D_cState Euler3D_pState::F(void) const {
   return (Euler3D_cState(d()*vx(), d()*sqr(vx()) + p(), 
                          d()*vx()*vy(), d()*vx()*vz(),
                          vx()*H()));
}

inline Euler3D_cState Euler3D_pState::F(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d()*W.vx(), W.d()*sqr(W.vx()) + W.p(),
                          W.d()*W.vx()*W.vy(), W.d()*W.vx()*W.vz(),
                          W.vx()*W.H()));
}

inline Euler3D_cState F(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d()*W.vx(), W.d()*sqr(W.vx()) + W.p(),
                          W.d()*W.vx()*W.vy(), W.d()*W.vx()*W.vz(),
                          W.vx()*W.H()));
}

inline Euler3D_cState Euler3D_pState::Fx(void) {
   return (Euler3D_cState(d()*vx(), d()*sqr(vx()) + p(), 
                          d()*vx()*vy(), d()*vx()*vz(), 
                          vx()*H()));
}

inline Euler3D_cState Euler3D_pState::Fx(void) const {
   return (Euler3D_cState(d()*vx(), d()*sqr(vx()) + p(),
                          d()*vx()*vy(), d()*vx()*vz(), 
                          vx()*H()));
}

inline Euler3D_cState Euler3D_pState::Fx(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d()*W.vx(), W.d()*sqr(W.vx()) + W.p(),
                          W.d()*W.vx()*W.vy(), W.d()*W.vx()*W.vz(), 
                          W.vx()*W.H()));
}

inline Euler3D_cState Fx(const Euler3D_pState &W) {
   return (Euler3D_cState(W.d()*W.vx(), W.d()*sqr(W.vx()) + W.p(),
                          W.d()*W.vx()*W.vy(), W.d()*W.vx()*W.vz(), 
                          W.vx()*W.H()));
}


/********************************************************
 * Useful 3D Euler state constants.                     *
 ********************************************************/
const Euler3D_pState Euler3D_W_STDATM(DENSITY_STDATM,
				      Vector3D_ZERO, PRESSURE_STDATM);
const Euler3D_pState Euler3D_W_VACUUM(ZERO, Vector3D_ZERO, ZERO);
const Euler3D_cState Euler3D_U_STDATM(Euler3D_W_STDATM);
const Euler3D_cState Euler3D_U_VACUUM(Euler3D_W_VACUUM);

extern Euler3D_pState RoeAverage(const Euler3D_pState &Wl,
	      	                 const Euler3D_pState &Wr);
extern Euler3D_cState FluxHLLE_x(const Euler3D_pState &Wl,
	      	                 const Euler3D_pState &Wr);
extern Euler3D_cState FluxHLLE_x(const Euler3D_cState &Ul,
	      	                 const Euler3D_cState &Ur);
extern Euler3D_cState FluxHLLE_n(const Euler3D_pState &Wl,
	      	                 const Euler3D_pState &Wr,
                                 const Vector3D &norm_dir);
extern Euler3D_cState FluxHLLE_n(const Euler3D_cState &Ul,
	      	                 const Euler3D_cState &Ur,
                                 const Vector3D &norm_dir);
extern Euler3D_pState Reflect(const Euler3D_pState &W,
	      	              const Vector3D &norm_dir);

#endif /* _EULER3D_STATE_INCLUDED  */

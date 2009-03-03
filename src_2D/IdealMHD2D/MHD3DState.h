/*!\file MHD3DState.h
  \brief Header file defining 3D/2D Ideal MHD Solution State Classes. */

#ifndef _MHD3D_STATE_INCLUDED
#define _MHD3D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"  /* Include math macro header files. */
#include "../CFD/CFD.h"    /* Include CFD header files. */
#include "../Math/Matrix.h"  /* Include matrix header files. */
#include "../Math/Vector2D.h" /* Include vector 2D header files. */
#include "../Math/Vector3D.h" /* Include vector 3D header files. */
#include "../Physics/GasConstants.h" /* Include gas constant header files. */
#include "../Physics/SolidConstants.h" /* Include solid constant header files. */

/* Define the classes. */

#define	NUM_VAR_MHD3D    8

class MHD3D_cState;

/*!
 * \class MHD3D_pState
 *
 * @brief Primitive solution state class definition for the 3D ideal MHD equations.
 *
 * Primitive solution state class definition for the 3D ideal MHD equations.
 * A 2D MHD form can be obtained by setting vz and Bz components to zero.
 *
 * \verbatim
 *                                                      
 * Member functions                                     
 *      d       -- Return density.                      
 *      v       -- Return flow velocity.                
 *      p       -- Return pressure.                     
 *      B1      -- Return perturbative magnetic field.  
 *      B0      -- Return intrinsic magnetic field.     
 *      g       -- Return specific heat ratio.          
 *      gm1     -- Return g-1                           
 *      gm1i    -- Return 1/(g-1).                      
 *      setgas  -- Set gas constants.                   
 *      B       -- Return total magnetic field.         
 *      T       -- Return temperature.                  
 *      e       -- Return specific internal energy.     
 *      E       -- Return total energy.                 
 *      E1      -- Return total perturbative energy.    
 *      h       -- Return specific total enthalpy.            
 *      h1      -- Return spec. perturbative enthalpy.  
 *      H       -- Return total enthalpy.               
 *      H1      -- Return total perturbative enthalpy.  
 *      a       -- Return sound speed.                  
 *      a2      -- Return sound speed square.           
 *      Va      -- Return Alfven wave velocity.         
 *      Va2     -- Return Alfven wave speed square.     
 *      s       -- Return specific entropy.             
 *      dv      -- Return momentum.                     
 *      U       -- Return conserved solution state.     
 *      F       -- Return solution flux.                
 *      lambda  -- Return eigenvalue.                   
 *      rp      -- Return primitive right eigenvector.  
 *      rc      -- Return conserved right eigenvector.  
 *      lp      -- Return primitive left eigenvector.   
 *      lc      -- Return conserved left eigenvector.    
 *                                                      
 * Member operators                                     
 *      W -- a primitive solution state                 
 *      c -- a scalar (double)                          
 *                                                      
 * W = W;                                               
 * c = W[i];                                            
 * W = W + W;                                           
 * W = W - W;                                           
 * c = W * W; (inner product)                           
 * W = c * W;                                           
 * W = W * c;                                           
 * W = W / c;                                           
 * W = +W;                                              
 * W = -W;                                              
 * W += W;                                              
 * W -= W;                                              
 * W == W;                                              
 * W != W;                                              
 * cout << W; (output function)                         
 * cin  >> W; (input function)                          
 * \endverbatim                                                     
 */
class MHD3D_pState{
public:

  //! @name Defined datatypes
  //@{
  // NONE
  //@}

  //@{ @name Associated constants to primitive variables:
  static double      g; //!< Specific heat ratio.
  static double    gm1; //!< g-1
  static double   gm1i; //!< 1/(g-1)
  //@}

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  MHD3D_pState(void):density(ONE), velocity(ZERO), pressure(ONE), PerturbativeB(ZERO), IntrinsicB(ZERO){ };

  //! Value Constructor
  explicit MHD3D_pState(const double &Val):density(Val), velocity(Val), pressure(Val), PerturbativeB(Val), IntrinsicB(Val){ };

  //! Copy constructor.
  MHD3D_pState(const MHD3D_pState &W):density(W.density), velocity(W.velocity), pressure(W.pressure),
				      PerturbativeB(W.PerturbativeB), IntrinsicB(W.IntrinsicB){ };

  //! Constructor with the conserved state
  MHD3D_pState(const MHD3D_cState &U);

  //! Constructor 2D with density, velocity, pressure and perturbative magnetic field 
  MHD3D_pState(const double &rho, const Vector2D &V_,
	       const Vector2D &Bm1, const double &pre) : density(rho), velocity(V_), pressure(pre), 
							 PerturbativeB(Bm1), IntrinsicB(ZERO) { };

  //! Constructor 3D with density, velocity, pressure and perturbative magnetic field 
  MHD3D_pState(const double &rho, const Vector3D &V_,
	       const Vector3D &Bm1, const double &pre) : density(rho), velocity(V_), pressure(pre), 
							 PerturbativeB(Bm1), IntrinsicB(ZERO) { };

  //! Constructor 2D with density, velocity, pressure, perturbative and intrinsic magnetic fields
  MHD3D_pState(const double &rho, const Vector2D &V_,
	       const Vector2D &Bm1, const Vector2D &Bm0,
	       const double &pre): density(rho), velocity(V_), pressure(pre),
				   PerturbativeB(Bm1), IntrinsicB(Bm0) { };

  //! Constructor 3D with density, velocity, pressure, perturbative and intrinsic magnetic fields
  MHD3D_pState(const double &rho, const Vector3D &V_,
	       const Vector3D &Bm1, const Vector3D &Bm0,
	       const double &pre): density(rho), velocity(V_), pressure(pre),
				   PerturbativeB(Bm1), IntrinsicB(Bm0) { };

  //! Constructor 2D with all components set and zero intrinsic magnetic field
  MHD3D_pState(const double &rho,
	       const double &vx,  const double &vy,
	       const double &B1x, const double &B1y,
	       const double &pre): density(rho), velocity(vx,vy), pressure(pre),
				   PerturbativeB(B1x,B1y), IntrinsicB(ZERO) { };

  //! Constructor 2D with all components set
  MHD3D_pState(const double &rho,
	       const double &vx,  const double &vy,
	       const double &B1x, const double &B1y,
	       const double &B0x, const double &B0y,
	       const double &pre) : density(rho), velocity(vx,vy), pressure(pre),
				    PerturbativeB(B1x,B1y), IntrinsicB(B0x,B0y){ };

  //! Constructor 3D with all components set
  MHD3D_pState(const double &rho,
	       const double &vx,  const double &vy, const double &vz,
	       const double &B1x, const double &B1y, const double &B1z,
	       const double &B0x, const double &B0y, const double &B0z,
	       const double &pre) : density(rho), velocity(vx,vy,vz), pressure(pre),
				    PerturbativeB(B1x,B1y,B1z), IntrinsicB(B0x,B0y,B0z){ };
  //@}

  //! @name Useful constant states
  //@{
  static const MHD3D_pState MHD3D_W_REF;
  static const MHD3D_pState MHD3D_W_ZERO;
  //@}

  //! @name Field access
  //@{
  double& d(void) {return density; }                     //!< Return density variable
  const double& d(void) const {return density; }         //!< Return density (read-only)

  Vector3D& v(void) {return velocity; }			 //!< Return velocity vector
  const Vector3D& v(void) const {return velocity; }	 //!< Return velocity vector (read-only)
  double& vx(void) {return velocity.x; }	         //!< Return x-direction velocity
  const double& vx(void) const {return velocity.x; }	 //!< Return x-direction velocity (read-only)
  double& vy(void) {return velocity.y; }		 //!< Return y-direction velocity
  const double& vy(void) const {return velocity.y; }	 //!< Return y-direction velocity (read-only)
  double& vz(void) {return velocity.z; }		 //!< Return z-direction velocity
  const double& vz(void) const {return velocity.z; }	 //!< Return z-direction velocity (read-only)

  double& p(void) {return pressure; }                    //!< Return pressure
  const double& p(void) const {return pressure; }        //!< Return pressure (read-only)

  Vector3D& B1(void) {return PerturbativeB; }		 //!< Return perturbative magnetic field
  const Vector3D& B1(void) const {return PerturbativeB;} //!< Return perturbative magnetic field (read-only)
  double& B1x(void) {return PerturbativeB.x; }		 //!< Return x-direction perturbative magnetic field
  const double& B1x(void) const {return PerturbativeB.x;}//!< Return x-direction perturbative magnetic field (read-only)
  double& B1y(void) {return PerturbativeB.y; }		 //!< Return y-direction perturbative magnetic field
  const double& B1y(void) const {return PerturbativeB.y;}//!< Return y-direction perturbative magnetic field (read-only)
  double& B1z(void) {return PerturbativeB.z; }		 //!< Return z-direction perturbative magnetic field
  const double& B1z(void) const {return PerturbativeB.z;}//!< Return z-direction perturbative magnetic field (read-only)

  Vector3D& B0(void) {return IntrinsicB; }		 //!< Return intrinsic magnetic field
  const Vector3D& B0(void) const {return IntrinsicB; }	 //!< Return intrinsic magnetic field (read-only)
  double& B0x(void) {return IntrinsicB.x; }		 //!< Return x-direction intrinsic magnetic field
  const double& B0x(void) const {return IntrinsicB.x;}   //!< Return x-direction intrinsic magnetic field (read-only)
  double& B0y(void) {return IntrinsicB.y; }		 //!< Return y-direction intrinsic magnetic field
  const double& B0y(void) const {return IntrinsicB.y;}   //!< Return y-direction intrinsic magnetic field (read-only)
  double& B0z(void) {return IntrinsicB.z; }		 //!< Return z-direction intrinsic magnetic field
  const double& B0z(void) const {return IntrinsicB.z;}   //!< Return z-direction intrinsic magnetic field (read-only)
  //@}

  //! Set gas constants.
  void setgas(void);
  void setgas(char *string_ptr);

  //! @name Physical quantities derived from the primitive variables
  //@{
  Vector3D B(void) const;                    //!< Return total magnetic field
  double T(void) const;                      //!< Return temperature
  double e(void) const;                      //!< Return specific internal energy
  double E(void) const;                      //!< Return total energy
  double E1(void) const;                     //!< Return total perturbative energy
  double h(void) const;	                     //!< Return specific total enthalpy
  double h1(void) const;                     //!< Return specific perturbative enthalpy
  double H(void) const;	                     //!< Return total enthalpy
  double H1(void) const;                     //!< Return total perturbative enthalpy
  double a(void) const;	                     //!< Return sound speed
  double a2(void) const;                     //!< Return sound speed squared
  Vector3D Va(void) const;                   //!< Return Alfven wave velocity
  double Va2(void) const;                    //!< Return Alfven wave speed squared
  double s(void) const;	                     //!< Return specific entropy
  Vector3D dv(void) const;	             //!< Return momentum vector
  double dv(const Vector2D &n) const;        //!< Return momentum projection along a given planar direction
  double dv(const Vector3D &n) const;        //!< Return momentum projection along a given 3D direction
  //@}

  MHD3D_cState U(void) const;		        //!< Return corresponding conserved solution state
  friend MHD3D_cState U(const MHD3D_pState &W); //!< Return the conserved solution state for a given primitive one
    
  MHD3D_cState F(void) const;	                //!< Return solution flux
  friend MHD3D_cState F(const MHD3D_pState &W); //!< Return solution flux for a given primitive state

  //! @name Eigenvalue(s) and Eigenvector(s)
  //@{
  MHD3D_pState lambda(void) const;      //!< Return the eigenvalues
  friend MHD3D_pState lambda(const MHD3D_pState &W){ return W.lambda(); }
  double lambda(int index) const;       //!< Return a particular eigenvalue
  friend double lambda(const MHD3D_pState &W, int index){ return W.lambda(index); }

  MHD3D_pState rp(int index) const;     //!< Primitive right eigenvector
  friend MHD3D_pState rp(const MHD3D_pState &W, int index){ return W.rp(index); }

  MHD3D_cState rc(int index) const;     //!< Conserved right eigenvector
  friend MHD3D_cState rc(const MHD3D_pState &W, int index);

  MHD3D_pState lp(int index) const;     //!< Primitive left eigenvector.
  friend MHD3D_pState lp(const MHD3D_pState &W, int index){ return W.lp(index); }

  MHD3D_cState lc(int index) const;     //!< Conserved left eigenvector 
  friend MHD3D_cState lc(const MHD3D_pState &W, int index);
  //@}

  //! Index operator (i.e. assigns an index to each variable)
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_MHD3D );
    switch(index) {
    case 1 :
      return density;
    case 2 :
      return velocity.x;
    case 3 :
      return velocity.y;
    case 4 :
      return velocity.z;
    case 5 :
      return PerturbativeB.x;
    case 6 :
      return PerturbativeB.y;
    case 7 :
      return PerturbativeB.z;
    case 8 :
      return pressure;
    default:
      return density;
    }
  }
    
  //! Read-only index operator
  const double &operator[](int index) const {
    assert( index >= 1 && index <= NUM_VAR_MHD3D );
    switch(index) {
    case 1 :
      return density;
    case 2 :
      return velocity.x;
    case 3 :
      return velocity.y;
    case 4 :
      return velocity.z;
    case 5 :
      return PerturbativeB.x;
    case 6 :
      return PerturbativeB.y;
    case 7 :
      return PerturbativeB.z;
    case 8 :
      return pressure;
    default:
      return density;
    }
  }

private:

  //! @name Solution state primitive variables:
  //@{
  double          density;       //!< Gas density.
  Vector3D        velocity;      //!< Flow velocity (3D vector).
  double          pressure;      //!< Gas pressure.
  Vector3D        PerturbativeB; //!< Perturbative magnetic field (3D vector).
  Vector3D        IntrinsicB;    //!< Intrinsic magnetic field (3D vector).
  //@}
  
};



/*!
 * Assign monatomic gas constants. 
 */
inline void MHD3D_pState::setgas(void) {
  g = GAMMA_MONATOMIC;
  gm1 = g - ONE;
  gm1i = ONE/gm1;
}

/*!
 * Assign gas constants. 
 */
inline void MHD3D_pState::setgas(char *string_ptr) {
  if (strcmp(string_ptr, "MONATOMIC") == 0) {
    g = GAMMA_MONATOMIC;
  } else if (strcmp(string_ptr, "DIATOMIC") == 0) {
    g = GAMMA_DIATOMIC;
  } else if (strcmp(string_ptr, "POLYATOMIC") == 0) {
    g = GAMMA_POLYATOMIC;
  } else if (strcmp(string_ptr, "BRIOWU") == 0) {
    g = GAMMA_TWO;
  } else if (strcmp(string_ptr, "ISOTHERMAL") == 0) {
    g = GAMMA_ISOTHERMAL;
  } /* endif */
  gm1 = g - ONE;
  gm1i = ONE/gm1;
}

/*!
 * Get total magnetic field. 
 */
inline Vector3D MHD3D_pState::B(void) const {
  return (B1()+B0());
}

/*!
 * Get temperature.
 */
inline double MHD3D_pState::T(void) const {
  return (pressure/density);
}

/*!
 * Specific internal energy.
 */
inline double MHD3D_pState::e(void) const {
  return (pressure/(gm1*density));
}

/*!
 * Total energy.
 */
inline double MHD3D_pState::E(void) const {
  return (pressure*gm1i + HALF*density*velocity.sqr() + HALF*B().sqr());
}

/*!
 * Total perturbative energy.
 */
inline double MHD3D_pState::E1(void) const {
  return (p()*gm1i + HALF*d()*v().sqr() + HALF*B1().sqr());
}

/*!
 * Specific total enthalpy. (i.e. it includes kinetic energy)
 */
inline double MHD3D_pState::h(void) const {
  return (g*p()/(gm1*d()) + HALF*v().sqr() + B().sqr()/d());
}

/*!
 * Specific total perturbative enthalpy.
 */
inline double MHD3D_pState::h1(void) const {
  return (g*p()/(gm1*d()) + HALF*v().sqr() + B1().sqr()/d());
}

/*!
 * Total enthalpy.
 */
inline double MHD3D_pState::H(void) const {
  return (g*gm1i*p() + HALF*d()*v().sqr() + B().sqr());
}

/*!
 * Total perturbative enthalpy.
 */
inline double MHD3D_pState::H1(void) const {
  return (g*gm1i*p() + HALF*d()*v().sqr() + B1().sqr());
}

/*!
 * Sound speed.
 */
inline double MHD3D_pState::a(void) const {
  return (sqrt(g*pressure/density));
}

/*!
 * Sound speed squared.
 */
inline double MHD3D_pState::a2(void) const {
  return (g*pressure/density);
}

/*!
 * Alfven wave velocity.
 */
inline Vector3D MHD3D_pState::Va(void) const {
  return (B()/sqrt(density));
}

/*!
 * Alfven wave speed squared.
 */
inline double MHD3D_pState::Va2(void) const {
  return ((B().sqr())/density);
}

/*!
 * Specific entropy.
 */
inline double MHD3D_pState::s(void) const {
  return (-gm1i*log(pressure/pow(density, g)));
}

/*!
 * Momentum vector
 */
inline Vector3D MHD3D_pState::dv(void) const {
  return density*velocity;
}

/*!
 * Momentum for a given plannar direction
 */
inline double MHD3D_pState::dv(const Vector2D & n) const {
  return density*(velocity*Vector3D(n));
}

/*!
 * Momentum for a given 3D direction
 */
inline double MHD3D_pState::dv(const Vector3D & n) const {
  return density*(velocity*n);
}

/*!
 * Determine eigenvalue(s) for a one-dimensional problem.
 */
inline MHD3D_pState MHD3D_pState::lambda(void) const {
  double c = a(), c2 = a2(), v1, v2, cs, cf;
  Vector3D ca = Va();
  v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
  cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
  cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
  return (MHD3D_pState(v().x - cf, v().x - ca.x, v().x - cs, v().x, v().x,
		       v().x + cs, v().x + ca.x, v().x + cf));
}

/*!
 * Associate an index to each eigenvalue for a one-dimensional problem.
 */
inline double MHD3D_pState::lambda(int index) const {
  double c, c2, v1, v2, cs, cf;
  assert( index >= 1 && index <= NUM_VAR_MHD3D );
  switch(index) {
  case 1 :
    c = a(); c2 = a2();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    return (v().x-cf);
  case 2 :
    return (v().x-Va().x);
  case 3 :
    c = a(); c2 = a2();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    return (v().x-cs);
  case 4 :
    return (v().x);
  case 5 :
    return (v().x);
  case 6 :
    c = a(); c2 = a2();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    return (v().x+cs);
  case 7 :
    return (v().x+Va().x);
  case 8 :
    c = a(); c2 = a2();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    return (v().x+cf);
  default:
    return (v().x);
  }
}

/*!
 * Primitive right eigenvector corresponding to each eigenvalue.
 */
inline MHD3D_pState MHD3D_pState::rp(int index) const {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD3D );
  switch(index) {
  case 1 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(alpha_f*d(), -alpha_f*cf, alpha_s*cs*beta_y*sgn(ca.x),
			 alpha_s*cs*beta_z*sgn(ca.x), ZERO,
			 alpha_s*c*beta_y*sqrt(d()), alpha_s*c*beta_z*sqrt(d()),
			 alpha_f*g*p()));
  case 2 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			 -beta_z*sqrt(d()/TWO), beta_y*sqrt(d()/TWO), ZERO));
  case 3 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(alpha_s*d(), -alpha_s*cs, -alpha_f*cf*beta_y*sgn(ca.x),
			 -alpha_f*cf*beta_z*sgn(ca.x), ZERO,
			 -alpha_f*c*beta_y*sqrt(d()), -alpha_f*c*beta_z*sqrt(d()),
			 alpha_s*g*p()));
  case 4 :
    return (MHD3D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  case 5 :
    return (MHD3D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
  case 6 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(alpha_s*d(), alpha_s*cs, alpha_f*cf*beta_y*sgn(ca.x),
			 alpha_f*cf*beta_z*sgn(ca.x), ZERO,
			 -alpha_f*c*beta_y*sqrt(d()), -alpha_f*c*beta_z*sqrt(d()),
			 alpha_s*g*p()));
  case 7 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			 beta_z*sqrt(d()/TWO), -beta_y*sqrt(d()/TWO), ZERO));
  case 8 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(alpha_f*d(), alpha_f*cf, -alpha_s*cs*beta_y*sgn(ca.x),
			 -alpha_s*cs*beta_z*sgn(ca.x), ZERO,
			 alpha_s*c*beta_y*sqrt(d()), alpha_s*c*beta_z*sqrt(d()),
			 alpha_f*g*p()));
  default: 
    return (MHD3D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  }
}

/*!
 * Primitive left eigenvector.
 */
inline MHD3D_pState MHD3D_pState::lp(int index) const {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD3D );
  switch(index) {
  case 1 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, -alpha_f*cf/(TWO*c2),
			 alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
			 alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
			 alpha_s*beta_y/(TWO*c*sqrt(d())),
			 alpha_s*beta_z/(TWO*c*sqrt(d())), alpha_f/(TWO*d()*c2)));
  case 2 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			 -beta_z/sqrt(TWO*d()), beta_y/sqrt(TWO*d()), ZERO));
  case 3 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, -alpha_s*cs/(TWO*c2),
			 -alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
			 -alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
			 -alpha_f*beta_y/(TWO*c*sqrt(d())),
			 -alpha_f*beta_z/(TWO*c*sqrt(d())), alpha_s/(TWO*d()*c2)));
  case 4 :
    return (MHD3D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, -ONE/a2()));
  case 5 :
    return (MHD3D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
  case 6 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, alpha_s*cs/(TWO*c2),
			 alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
			 alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
			 -alpha_f*beta_y/(TWO*c*sqrt(d())),
			 -alpha_f*beta_z/(TWO*c*sqrt(d())), alpha_s/(TWO*d()*c2)));
  case 7 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			 beta_z/sqrt(TWO*d()), -beta_y/sqrt(TWO*d()), ZERO));
  case 8 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_pState(ZERO, alpha_f*cf/(TWO*c2),
			 -alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
			 -alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
			 alpha_s*beta_y/(TWO*c*sqrt(d())),
			 alpha_s*beta_z/(TWO*c*sqrt(d())), alpha_f/(TWO*d()*c2)));
  default: 
    return (MHD3D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  }
}



















/*!
 * \class MHD3D_cState                                  
 *                        
 * @brief Conserved solution state class definition for the 3D ideal MHD equations.
 *                              
 * Member functions                                     
 *      d       -- Return density.                      
 *      dv      -- Return momentum.                     
 *      B1      -- Return perturbative magentic field.  
 *      B0      -- Return intrinsic magnetic field.     
 *      E1      -- Return total perturbative energy.    
 *      g       -- Return specific heat ratio.          
 *      gm1     -- Return g-1.                          
 *      gm1i    -- Return 1/(g-1).                      
 *      setgas  -- Set gas constants.                   
 *      v       -- Return flow velocity.                
 *      p       -- Return pressure.                     
 *      B       -- Return total magnetic field.         
 *      T       -- Return temperature.                  
 *      e       -- Return specific internal energy.     
 *      E       -- Return total energy.                 
 *      h       -- Return specific enthalpy.            
 *      h1      -- Return perturbative enthalpy.        
 *      H       -- Return total enthalpy.               
 *      H1      -- Return total perturbative enthalpy.  
 *      a       -- Return sound speed.                  
 *      a2      -- Return sound speed square.           
 *      Va      -- Return Alfven wave velocity.         
 *      Va2     -- Return Alfven wave speed square.     
 *      s       -- Return specific entropy.             
 *      W       -- Return primitive solution state.     
 *                                                      
 * Member operators                                     
 *      U -- a primitive solution state                 
 *      c -- a scalar (double)                          
 *                                                      
 * U = U;                                               
 * c = U[i];                                            
 * U = U + U;                                           
 * U = U - U;                                           
 * c = U * U; (inner product)                           
 * U = c * U;                                           
 * U = U * c;                                           
 * U = U / c;                                           
 * U = +U;                                              
 * U = -U;                                              
 * U += U;                                              
 * U -= U;                                              
 * U == U;                                              
 * U != U;                                              
 * cout << U; (output function)                         
 * cin  >> U; (input function)                          
 *                                                      
 */
class MHD3D_cState{
public:

  //! @name Defined datatypes
  //@{
  // NONE
  //@}

  //@{ @name Associated constants to primitive variables:
  static double      g; //!< Specific heat ratio.
  static double    gm1; //!< g-1
  static double   gm1i; //!< 1/(g-1)
  //@}

  //! @name Useful constant states
  //@{
  static const MHD3D_cState MHD3D_U_REF;
  static const MHD3D_cState MHD3D_U_ZERO;
  //@}

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  MHD3D_cState(void): density(ONE), momentum(ZERO), PerturbativeB(ZERO), IntrinsicB(ZERO),
		      PerturbativeEnergy(ONE/(GAMMA_MONATOMIC-ONE)){ };

  //! Value Constructor
  explicit MHD3D_cState(const double &Val): density(Val), momentum(Val), PerturbativeB(Val), IntrinsicB(Val),
					    PerturbativeEnergy(Val){ };

  //! Copy constructor.
  MHD3D_cState(const MHD3D_cState &U): density(U.density), momentum(U.momentum),
				       PerturbativeB(U.PerturbativeB), IntrinsicB(U.IntrinsicB),
				       PerturbativeEnergy(U.PerturbativeEnergy){ };

  //! Constructor with the primitive state
  MHD3D_cState(const MHD3D_pState &W);

  //! Constructor 2D with density, momentum, perturbative magnetic field and total perturbative energy
  MHD3D_cState(const double &rho,
	       const Vector2D &rhoV,
	       const Vector2D &Bm1,
	       const double &Etotal): density(rho), momentum(rhoV),
				      PerturbativeB(Bm1), IntrinsicB(ZERO), PerturbativeEnergy(Etotal){ };

  //! Constructor 3D with density, momentum, perturbative magnetic field and total perturbative energy
  MHD3D_cState(const double &rho,
	       const Vector3D &rhoV,
	       const Vector3D &Bm1,
	       const double &Etotal): density(rho), momentum(rhoV),
				      PerturbativeB(Bm1), IntrinsicB(ZERO), PerturbativeEnergy(Etotal){ };

  //! Constructor 2D with density, momentum, perturbative and intrinsic magnetic fields, total perturbative energy
  MHD3D_cState(const double &rho,
	       const Vector2D &rhoV,
	       const Vector2D &Bm1,
	       const Vector2D &Bm0,
	       const double &Etotal): density(rho), momentum(rhoV),
				      PerturbativeB(Bm1), IntrinsicB(Bm0), PerturbativeEnergy(Etotal){ };

  //! Constructor 3D with density, momentum, perturbative and intrinsic magnetic fields, total perturbative energy
  MHD3D_cState(const double &rho,
	       const Vector3D &rhoV,
	       const Vector3D &Bm1,
	       const Vector3D &Bm0,
	       const double &Etotal): density(rho), momentum(rhoV),
				      PerturbativeB(Bm1), IntrinsicB(Bm0), PerturbativeEnergy(Etotal){ };

  //! Constructor 2D with all components set and zero intrinsic magnetic field
  MHD3D_cState(const double &rho,
	       const double &rhovx,
	       const double &rhovy,
	       const double &_B1x_,
	       const double &_B1y_,
	       const double &Etotal): density(rho), momentum(rhovx,rhovy), 
				      PerturbativeB(_B1x_,_B1y_), IntrinsicB(ZERO), PerturbativeEnergy(Etotal){ };

  //! Constructor 2D with all components set
  MHD3D_cState(const double &rho,
	       const double &rhovx,
	       const double &rhovy,
	       const double &_B1x_,
	       const double &_B1y_,
	       const double &_B0x_,
	       const double &_B0y_,
	       const double &Etotal): density(rho), momentum(rhovx,rhovy),
				      PerturbativeB(_B1x_,_B1y_), IntrinsicB(_B0x_,_B0y_), PerturbativeEnergy(Etotal){ };

  //! Constructor 3D with all components set
  MHD3D_cState(const double &rho,
	       const double &rhovx,
	       const double &rhovy,
	       const double &rhovz,
	       const double &_B1x_,
	       const double &_B1y_,
	       const double &_B1z_,
	       const double &_B0x_,
	       const double &_B0y_,
	       const double &_B0z_,
	       const double &Etotal): density(rho), momentum(rhovx,rhovy,rhovz),
				      PerturbativeB(_B1x_,_B1y_,_B1z_), IntrinsicB(_B0x_,_B0y_,_B0z_),
				      PerturbativeEnergy(Etotal){ };
  //@}    


  //! @name Field access
  //@{
  double& d(void) {return density; }                     //!< Return density variable
  const double& d(void) const {return density; }         //!< Return density (read-only)

  Vector3D& dv(void) {return momentum; }		 //!< Return momentum vector
  const Vector3D& dv(void) const {return momentum; }	 //!< Return momentum vector (read-only)
  double& dvx(void) {return momentum.x; }	         //!< Return x-direction momentum
  const double& dvx(void) const {return momentum.x; }	 //!< Return x-direction momentum (read-only)
  double& dvy(void) {return momentum.y; }		 //!< Return y-direction momentum
  const double& dvy(void) const {return momentum.y; }	 //!< Return y-direction momentum (read-only)
  double& dvz(void) {return momentum.z; }		 //!< Return z-direction momentum
  const double& dvz(void) const {return momentum.z; }	 //!< Return z-direction momentum (read-only)
  
  Vector3D& B1(void) {return PerturbativeB; }		 //!< Return perturbative magnetic field
  const Vector3D& B1(void) const {return PerturbativeB;} //!< Return perturbative magnetic field (read-only)
  double& B1x(void) {return PerturbativeB.x; }		 //!< Return x-direction perturbative magnetic field
  const double& B1x(void) const {return PerturbativeB.x;}//!< Return x-direction perturbative magnetic field (read-only)
  double& B1y(void) {return PerturbativeB.y; }		 //!< Return y-direction perturbative magnetic field
  const double& B1y(void) const {return PerturbativeB.y;}//!< Return y-direction perturbative magnetic field (read-only)
  double& B1z(void) {return PerturbativeB.z; }		 //!< Return z-direction perturbative magnetic field
  const double& B1z(void) const {return PerturbativeB.z;}//!< Return z-direction perturbative magnetic field (read-only)

  Vector3D& B0(void) {return IntrinsicB; }		 //!< Return intrinsic magnetic field
  const Vector3D& B0(void) const {return IntrinsicB; }	 //!< Return intrinsic magnetic field (read-only)
  double& B0x(void) {return IntrinsicB.x; }		 //!< Return x-direction intrinsic magnetic field
  const double& B0x(void) const {return IntrinsicB.x;}   //!< Return x-direction intrinsic magnetic field (read-only)
  double& B0y(void) {return IntrinsicB.y; }		 //!< Return y-direction intrinsic magnetic field
  const double& B0y(void) const {return IntrinsicB.y;}   //!< Return y-direction intrinsic magnetic field (read-only)
  double& B0z(void) {return IntrinsicB.z; }		 //!< Return z-direction intrinsic magnetic field
  const double& B0z(void) const {return IntrinsicB.z;}   //!< Return z-direction intrinsic magnetic field (read-only)

  double& E1(void) {return PerturbativeEnergy; }              //!< Return total perturbative energy
  const double& E1(void) const {return PerturbativeEnergy; }  //!< Return total perturbative energy (read-only)
  //@}

  //! Set gas constants.
  void setgas(void);
  void setgas(char *string_ptr);

  //! @name Physical quantities derived from the conserved variables
  //@{
  Vector3D v(void) const;	         //!< Return flow velocity (3D vector).
  double v(const Vector2D &n) const;	 //!< Return flow velocity along a given planar direction.
  double v(const Vector3D &n) const;	 //!< Return flow velocity along a given 3D direction.
  double p(void) const;			 //!< Return gas pressure.
  Vector3D B(void) const;		 //!< Return total magnetic field. (3D vector)
  double T(void) const;			 //!< Return temparature.
  double e(void) const;			 //!< Return specific internal energy.
  double E(void) const;			 //!< Return total energy.
  double h(void) const;			 //!< Return specific enthalpy.
  double h1(void) const;		 //!< Return specific perturbative enthalpy.
  double H(void) const;			 //!< Return total enthalpy.
  double H1(void) const;		 //!< Return total perturbative enthalpy.
  double a(void) const;			 //!< Return sound speed.
  double a2(void) const;		 //!< Return sound speed squared.
  Vector3D Va(void) const;		 //!< Return Alfven wave velocity.
  double Va2(void) const;		 //!< Return Alfven wave speed squared.
  double s(void) const;			 //!< Return specific entropy.
  //@}

  MHD3D_pState W(void) const;	         //!< Return corresponding primitive solution state 
  friend MHD3D_pState W(const MHD3D_cState &U){ return U.W(); } //!< Return the primitive solution state for a given conserved one
    
  //! Index operator (i.e. assigns an index to each variable)
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_MHD3D );
    switch(index) {
    case 1 :
      return density;
    case 2 :
      return momentum.x;
    case 3 :
      return momentum.y;
    case 4 :
      return momentum.z;
    case 5 :
      return PerturbativeB.x;
    case 6 :
      return PerturbativeB.y;
    case 7 :
      return PerturbativeB.z;
    case 8 :
      return PerturbativeEnergy;
    default:
      return density;
    }
  }
    
  const double &operator[](int index) const {
    assert( index >= 1 && index <= NUM_VAR_MHD3D );
    switch(index) {
    case 1 :
      return density;
    case 2 :
      return momentum.x;
    case 3 :
      return momentum.y;
    case 4 :
      return momentum.z;
    case 5 :
      return PerturbativeB.x;
    case 6 :
      return PerturbativeB.y;
    case 7 :
      return PerturbativeB.z;
    case 8 :
      return PerturbativeEnergy;
    default:
      return density;
    }
  }


#if 0
  /* Binary arithmetic operators. */
  friend MHD3D_cState operator +(const MHD3D_cState &U1, const MHD3D_cState &U2);
  friend MHD3D_cState operator -(const MHD3D_cState &U1, const MHD3D_cState &U2);
  friend double operator *(const MHD3D_cState &U1, const MHD3D_cState &U2);
  friend MHD3D_cState operator *(const MHD3D_cState &U, const double &a);
  friend MHD3D_cState operator *(const double &a, const MHD3D_cState &U);
  friend MHD3D_cState operator /(const MHD3D_cState &U, const double &a);

  /* Unary arithmetic operators. */
  friend MHD3D_cState operator +(const MHD3D_cState &U);
  friend MHD3D_cState operator -(const MHD3D_cState &U);

  /* Shortcut arithmetic operators. */
  friend MHD3D_cState &operator +=(MHD3D_cState &U1, const MHD3D_cState &U2);
  friend MHD3D_cState &operator -=(MHD3D_cState &U1, const MHD3D_cState &U2);
    
  /* Relational operators. */
  friend int operator ==(const MHD3D_cState &U1, const MHD3D_cState &U2);
  friend int operator !=(const MHD3D_cState &U1, const MHD3D_cState &U2);
    
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file, const MHD3D_cState &U);
  friend istream &operator >> (istream &in_file,  MHD3D_cState &U);
#endif

private:
  
  //! @name Solution state conserved variables:
  //@{
  double         density;             //!< Gas density.
  Vector3D       momentum;            //!< Gas momentum (3D vector).
  Vector3D       PerturbativeB;       //!< Perturbative magnetic field (3D vector).
  Vector3D       IntrinsicB;          //!< Intrinsic magnetic field (3D vector).
  double         PerturbativeEnergy;  //!< Total perturbative energy.
  //@}

};


/*!
 * Constructor with the primitive state
 */
inline MHD3D_cState::MHD3D_cState(const MHD3D_pState &W){
  density = W.d();
  momentum = W.dv();
  PerturbativeB = W.B1();
  IntrinsicB = W.B0();
  PerturbativeEnergy = W.E1();
}

/*!
 * Assign monatomic gas constants. 
 */
inline void MHD3D_cState::setgas(void) {
  g = GAMMA_MONATOMIC;
  gm1 = g - ONE;
  gm1i = ONE/gm1;
}

/*!
 * Assign gas constants. 
 */
inline void MHD3D_cState::setgas(char *string_ptr) {
  if (strcmp(string_ptr, "MONATOMIC") == 0) {
    g = GAMMA_MONATOMIC;
  } else if (strcmp(string_ptr, "DIATOMIC") == 0) {
    g = GAMMA_DIATOMIC;
  } else if (strcmp(string_ptr, "POLYATOMIC") == 0) {
    g = GAMMA_POLYATOMIC;
  } else if (strcmp(string_ptr, "BRIOWU") == 0) {
    g = GAMMA_TWO;
  } else if (strcmp(string_ptr, "ISOTHERMAL") == 0) {
    g = GAMMA_ISOTHERMAL;
  } /* endif */
  gm1 = g - ONE;
  gm1i = ONE/gm1;
}

/*!
 * Get flow velocity.
 */
inline Vector3D MHD3D_cState::v(void) const {
  return (dv()/d());
}

/*!
 * Get flow velocity along a given direction in 2D.
 * Third coordinate is considered zero.
 */
inline double MHD3D_cState::v(const Vector2D &n) const {
  return ((dv()*Vector3D(n))/d());
}

/*!
 * Get flow velocity along a given direction in 3D
 */
inline double MHD3D_cState::v(const Vector3D &n) const {
  return ((dv()*n)/d());
}

/*!
 * Get gas pressure.
 */
inline double MHD3D_cState::p(void) const {
  return (gm1*(E1() - HALF*dv().sqr()/d() - HALF*B1().sqr()));
}

/*!
 * Get total magnetic field.
 */
inline Vector3D MHD3D_cState::B(void) const {
  return (B1()+B0());
}

/*!
 * Get gas temperature.
 */
inline double MHD3D_cState::T(void) const {
  return (p()/d());
}

/*!
 * Get specific internal energy.
 */
inline double MHD3D_cState::e(void) const {
  return ((E1() - HALF*dv().sqr()/d() - HALF*B1().sqr())/d());
}

/*!
 * Get total energy.
 */
inline double MHD3D_cState::E(void) const {
  return (E1() + HALF*B0().sqr() + B0()*B1());
}

/*!
 * Get specific enthalpy.
 */
inline double MHD3D_cState::h(void) const {
  return ( (g*E1() - gm1*HALF*dv().sqr()/d() - (g-TWO)*HALF*B1().sqr() + B0().sqr() + TWO*B0()*B1())/d() );
}

/*!
 * Get specific perturbative enthalpy.
 */
inline double MHD3D_cState::h1(void) const {
  return ( (g*E1() - gm1*HALF*dv().sqr()/d() - (g-TWO)*HALF*B1().sqr())/d() );
}

/*!
 * Get total enthalpy.
 */
inline double MHD3D_cState::H(void) const {
  return ( g*E1() - gm1*HALF*dv().sqr()/d() - (g-TWO)*HALF*B1().sqr() + B0().sqr() + TWO*B0()*B1() );
}

/*!
 * Get total perturbative enthalpy.
 */
inline double MHD3D_cState::H1(void) const {
  return ( g*E1() - gm1*HALF*dv().sqr()/d() - (g-TWO)*HALF*B1().sqr() );
}

/*!
 * Get sound speed.
 */
inline double MHD3D_cState::a(void) const {
  return ( sqrt(g*p()/d()) );
}

/*!
 * Get sound speed squared.
 */
inline double MHD3D_cState::a2(void) const {
  return ( g*p()/d() );
}

/*!
 * Get Alfven wave velocity.
 */
inline Vector3D MHD3D_cState::Va(void) const {
  return ( B()/sqrt(d()) );
}

/*!
 * Get Alfven wave speed squared.
 */
inline double MHD3D_cState::Va2(void) const {
  return ( B().sqr()/d() );
}

/*!
 * Get specific entropy.
 */
inline double MHD3D_cState::s(void) const {
  return (-gm1i*log(gm1*(E1() - HALF*dv().sqr()/d() - HALF*B1().sqr())/pow(d(), g)));
}



/*********************************************************************************
 *********************************************************************************
 *                                                                               *
 *             Member functions of MHD3D_pState that                             *
 *             need a completely defined MHD3D_cState                            *
 *                                                                               *
 *********************************************************************************
 ********************************************************************************/

/*! Constructor with the conserved state */
inline MHD3D_pState::MHD3D_pState(const MHD3D_cState &U){
  density = U.d();
  velocity = U.v();
  PerturbativeB = U.B1();
  IntrinsicB = U.B0();
  pressure = U.p();
}

/*! Return conserved solution state */
inline MHD3D_cState MHD3D_pState::U(void) const {
  return MHD3D_cState(density, dv(), B1(), B0(), E1());
}

/*! Return the conserved solution state for a given primitive one */
inline MHD3D_cState U(const MHD3D_pState &W){ return W.U(); }

/*! Solution flux */
inline MHD3D_cState MHD3D_pState::F(void) const {
  return (MHD3D_cState(d()*v().x,
		       d()*v().x*v().x - B1().x*B1().x - (B0().x*B1().x + B1().x*B0().x) + p() + HALF*B1().sqr() + B1()*B0(),
		       d()*v().x*v().y - B1().x*B1().y - (B0().x*B1().y + B1().x*B0().y),
		       d()*v().x*v().z - B1().x*B1().z - (B0().x*B1().z + B1().x*B0().z),
		       v().x*B1().x - B1().x*v().x + v().x*B0().x - B0().x*v().x,
		       v().x*B1().y - B1().x*v().y + v().x*B0().y - B0().x*v().y,
		       v().x*B1().z - B1().x*v().z + v().x*B0().z - B0().x*v().z,
		       v().x*H1() - (B1()*v())*B1().x + (B0()*B1())*v().x - (B1()*v())*B0().x));
}
    
/*! Return solution flux for a given primitive state */
inline MHD3D_cState F(const MHD3D_pState &W){ return W.F(); } 

/*! Conserved right eigenvector. */
inline MHD3D_cState MHD3D_pState::rc(int index) const {
  double c, c2, v1, v2, cs, cf,
    alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD3D );
  switch(index) {
  case 1 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = alpha_s*(c*sqrt(d())*(beta_y*B1().y+beta_z*B1().z)+
		     d()*cs*sgn(ca.x)*(v().y*beta_y+v().z*beta_z));
    return (MHD3D_cState(alpha_f*d(),
			 alpha_f*d()*(v().x-cf),
			 d()*(alpha_f*v().y+alpha_s*cs*beta_y*sgn(ca.x)),
			 d()*(alpha_f*v().z+alpha_s*cs*beta_z*sgn(ca.x)),
			 ZERO,
			 alpha_s*c*beta_y*sqrt(d()),
			 alpha_s*c*beta_z*sqrt(d()),
			 alpha_f*(HALF*d()*sqr(v())+g*p()*gm1i-d()*v().x*cf)+gamma));
  case 2 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_cState(ZERO,
			 ZERO,
			 -beta_z*d()/sqrt(TWO),
			 beta_y*d()/sqrt(TWO),
			 ZERO,
			 -beta_z*sqrt(d()/TWO),
			 beta_y*sqrt(d()/TWO),
			 ((v().z*beta_y-v().y*beta_z)-(B1().y*beta_z-B1().z*beta_y))*sqrt(d()/TWO)));
  case 3 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = alpha_f*(c*sqrt(d())*(beta_y*B1().y+beta_z*B1().z)+
		     d()*cf*sgn(ca.x)*(v().y*beta_y+v().z*beta_z));
    return (MHD3D_cState(alpha_s*d(),
			 alpha_s*d()*(v().x-cs),
			 d()*(alpha_s*v().y-alpha_f*cf*beta_y*sgn(ca.x)),
			 d()*(alpha_s*v().z-alpha_f*cf*beta_z*sgn(ca.x)),
			 ZERO,
			 -alpha_f*c*beta_y*sqrt(d()),
			 -alpha_f*c*beta_z*sqrt(d()),
			 alpha_s*(HALF*d()*sqr(v())+g*p()*gm1i-d()*v().x*cs)-gamma));
  case 4 :
    return (MHD3D_cState(ONE, v(), Vector3D_ZERO, HALF*sqr(v())));
  case 5 :
    return (MHD3D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, B1().x));
  case 6 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = alpha_f*(c*sqrt(d())*(beta_y*B1().y+beta_z*B1().z)-
		     d()*cf*sgn(ca.x)*(v().y*beta_y+v().z*beta_z));
    return (MHD3D_cState(alpha_s*d(),
			 alpha_s*d()*(v().x+cs),
			 d()*(alpha_s*v().y+alpha_f*cf*beta_y*sgn(ca.x)),
			 d()*(alpha_s*v().z+alpha_f*cf*beta_z*sgn(ca.x)),
			 ZERO,
			 -alpha_f*c*beta_y*sqrt(d()),
			 -alpha_f*c*beta_z*sqrt(d()),
			 alpha_s*(HALF*d()*sqr(v())+g*p()*gm1i+d()*v().x*cs)-gamma));
  case 7 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_cState(ZERO,
			 ZERO,
			 -beta_z*d()/sqrt(TWO),
			 beta_y*d()/sqrt(TWO),
			 ZERO,
			 beta_z*sqrt(d()/TWO),
			 -beta_y*sqrt(d()/TWO),
			 ((v().z*beta_y-v().y*beta_z)+(B1().y*beta_z-B1().z*beta_y))*sqrt(d()/TWO)));
  case 8 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = alpha_s*(c*sqrt(d())*(beta_y*B1().y+beta_z*B1().z)-
		     d()*cs*sgn(ca.x)*(v().y*beta_y+v().z*beta_z));
    return (MHD3D_cState(alpha_f*d(),
			 alpha_f*d()*(v().x+cf),
			 d()*(alpha_f*v().y-alpha_s*cs*beta_y*sgn(ca.x)),
			 d()*(alpha_f*v().z-alpha_s*cs*beta_z*sgn(ca.x)),
			 ZERO,
			 alpha_s*c*beta_y*sqrt(d()),
			 alpha_s*c*beta_z*sqrt(d()),
			 alpha_f*(HALF*d()*sqr(v())+g*p()*gm1i+d()*v().x*cf)+gamma));
  default: 
    return (MHD3D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  }
}

/*! Return the conserved right eigenvector */
inline MHD3D_cState rc(const MHD3D_pState &W, int index){ return W.rc(index); }

/*! Conserved left eigenvector. */
inline MHD3D_cState MHD3D_pState::lc(int index) const {
  double c, c2, v1, v2, cs, cf,
    alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD3D );
  switch(index) {
  case 1 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = -HALF*alpha_s*cs*sgn(ca.x)*(v().y*beta_y+v().z*beta_z)/(d()*c2);
    return (MHD3D_cState(HALF*alpha_f*(HALF*gm1*sqr(v())+v().x*cf)/(d()*c2)+gamma,
			 -HALF*alpha_f*(gm1*v().x+cf)/(d()*c2),
			 -HALF*(gm1*alpha_f*v().y-alpha_s*cs*beta_y*sgn(ca.x))/(d()*c2),
			 -HALF*(gm1*alpha_f*v().z-alpha_s*cs*beta_z*sgn(ca.x))/(d()*c2),
			 -HALF*gm1*alpha_f*B1().x/(d()*c2),
			 HALF*(alpha_s*c*beta_y*sqrt(d())-gm1*alpha_f*B1().y)/(d()*c2),
			 HALF*(alpha_s*c*beta_z*sqrt(d())-gm1*alpha_f*B1().z)/(d()*c2),
			 HALF*gm1*alpha_f/(d()*c2)));
  case 2 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_cState((v().y*beta_z - v().z*beta_y)/(sqrt(TWO)*d()),
			 ZERO,
			 -beta_z/(sqrt(TWO)*d()),
			 beta_y/(sqrt(TWO)*d()),
			 ZERO,
			 -beta_z/sqrt(TWO*d()),
			 beta_y/sqrt(TWO*d()),
			 ZERO));
  case 3 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = HALF*alpha_f*cf*sgn(ca.x)*(v().y*beta_y+v().z*beta_z)/(d()*c2);
    return (MHD3D_cState(HALF*alpha_s*(HALF*gm1*sqr(v())+v().x*cs)/(d()*c2)+gamma,
			 -HALF*alpha_s*(gm1*v().x+cs)/(d()*c2),
			 -HALF*(gm1*alpha_s*v().y+alpha_f*cf*beta_y*sgn(ca.x))/(d()*c2),
			 -HALF*(gm1*alpha_s*v().z+alpha_f*cf*beta_z*sgn(ca.x))/(d()*c2),
			 -HALF*gm1*alpha_s*B1().x/(d()*c2),
			 -HALF*(alpha_f*c*beta_y*sqrt(d())+gm1*alpha_s*B1().y)/(d()*c2),
			 -HALF*(alpha_f*c*beta_z*sqrt(d())+gm1*alpha_s*B1().z)/(d()*c2),
			 HALF*gm1*alpha_s/(d()*c2)));
  case 4 :
    c2 = a2();
    return (MHD3D_cState(ONE-HALF*gm1*sqr(v())/c2,
			 gm1*v().x/c2,
			 gm1*v().y/c2,
			 gm1*v().z/c2,
			 gm1*B1().x/c2,
			 gm1*B1().y/c2,
			 gm1*B1().z/c2,
			 -gm1/c2));
  case 5 :
    return (MHD3D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
  case 6 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = -HALF*alpha_f*cf*sgn(ca.x)*(v().y*beta_y+v().z*beta_z)/(d()*c2);
    return (MHD3D_cState(HALF*alpha_s*(HALF*gm1*sqr(v())-v().x*cs)/(d()*c2)+gamma,
			 -HALF*alpha_s*(gm1*v().x-cs)/(d()*c2),
			 -HALF*(gm1*alpha_s*v().y-alpha_f*cf*beta_y*sgn(ca.x))/(d()*c2),
			 -HALF*(gm1*alpha_s*v().z-alpha_f*cf*beta_z*sgn(ca.x))/(d()*c2),
			 -HALF*gm1*alpha_s*B1().x/(d()*c2),
			 -HALF*(alpha_f*c*beta_y*sqrt(d())+gm1*alpha_s*B1().y)/(d()*c2),
			 -HALF*(alpha_f*c*beta_z*sqrt(d())+gm1*alpha_s*B1().z)/(d()*c2),
			 HALF*gm1*alpha_s/(d()*c2)));
  case 7 :
    ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    return (MHD3D_cState((v().y*beta_z - v().z*beta_y)/(sqrt(TWO)*d()),
			 ZERO,
			 -beta_z/(sqrt(TWO)*d()),
			 beta_y/(sqrt(TWO)*d()),
			 ZERO,
			 beta_z/sqrt(TWO*d()),
			 -beta_y/sqrt(TWO*d()),
			 ZERO));
  case 8 :
    c = a(); c2 = a2(); ca = Va();
    v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
    cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
    cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
    v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
    if (v1 >= TOLER) {
      alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
    } else if (ca.x <= c) {
      alpha_f = ONE; alpha_s = ZERO;
    } else {
      alpha_f = ZERO; alpha_s = ONE;
    } /* endif */
    if (v2 >= TOLER) {
      beta_y = ca.y/v2; beta_z = ca.z/v2;
    } else {
      beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
    } /* endif */
    gamma = HALF*alpha_s*cs*sgn(ca.x)*(v().y*beta_y+v().z*beta_z)/(d()*c2);
    return (MHD3D_cState(HALF*alpha_f*(HALF*gm1*sqr(v())-v().x*cf)/(d()*c2)+gamma,
			 -HALF*alpha_f*(gm1*v().x-cf)/(d()*c2),
			 -HALF*(gm1*alpha_f*v().y+alpha_s*cs*beta_y*sgn(ca.x))/(d()*c2),
			 -HALF*(gm1*alpha_f*v().z+alpha_s*cs*beta_z*sgn(ca.x))/(d()*c2),
			 -HALF*gm1*alpha_f*B1().x/(d()*c2),
			 HALF*(alpha_s*c*beta_y*sqrt(d())-gm1*alpha_f*B1().y)/(d()*c2),
			 HALF*(alpha_s*c*beta_z*sqrt(d())-gm1*alpha_f*B1().z)/(d()*c2),
			 HALF*gm1*alpha_f/(d()*c2)));
  default: 
    return (MHD3D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  }
}

/*! Return the conserved left eigenvector */
inline MHD3D_cState lc(const MHD3D_pState &W, int index){ return W.lc(index); }


#endif  /* _MHD3D_STATE_INCLUDED */
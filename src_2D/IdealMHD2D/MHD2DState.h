/*!\file MHD2DState.h
  \brief Header file defining 2D Ideal MHD Solution State Classes. */

#ifndef _MHD2D_STATE_INCLUDED
#define _MHD2D_STATE_INCLUDED

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
#include "../Physics/GasConstants.h" /* Include gas constant header files. */
#include "../Physics/SolidConstants.h" /* Include solid constant header files. */

/* Define the classes. */

#define	NUM_VAR_MHD2D    6

class MHD2D_cState;

/*!
 * \class MHD2D_pState
 *
 * @brief Primitive solution state class definition for the 2D ideal MHD equations.
 *
 * Primitive solution state class definition for the 2D ideal MHD equations.
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
 *      h       -- Return specific enthalpy.            
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
class MHD2D_pState{
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
  MHD2D_pState(void):density(ONE), velocity(ZERO), pressure(ONE), PerturbativeB(ZERO), IntrinsicB(ZERO){ };

  //! Value Constructor
  explicit MHD2D_pState(const double &Val):density(Val), velocity(Val), pressure(Val), PerturbativeB(Val), IntrinsicB(Val){ };

  //! Copy constructor.
  MHD2D_pState(const MHD2D_pState &W):density(W.density), velocity(W.velocity), pressure(W.pressure),
				      PerturbativeB(W.PerturbativeB), IntrinsicB(W.IntrinsicB){ };

  //! Constructor with the conserved state
  MHD2D_pState(const MHD2D_cState &U);

  //! Constructor with density, velocity, pressure and perturbative magnetic field
  MHD2D_pState(const double &rho, const Vector2D &V_,
	       const Vector2D &Bm1, const double &pre) : density(rho), velocity(V_), pressure(pre), 
							 PerturbativeB(Bm1) { B0().zero(); }

  //! Constructor with density, velocity, pressure, perturbative and intrinsic magnetic fields
  MHD2D_pState(const double &rho, const Vector2D &V_,
	       const Vector2D &Bm1, const Vector2D &Bm0,
	       const double &pre): density(rho), velocity(V_), pressure(pre),
				   PerturbativeB(Bm1), IntrinsicB(Bm0) { };

  //! Constructor with all components set and zero intrinsic magnetic field
  MHD2D_pState(const double &rho,
	       const double &vx,  const double &vy,
	       const double &B1x, const double &B1y,
	       const double &pre): density(rho), velocity(vx,vy), pressure(pre),
				   PerturbativeB(B1x,B1y) { B0().zero(); }

  //! Constructor with all components set
  MHD2D_pState(const double &rho,
	       const double &vx,  const double &vy,
	       const double &B1x, const double &B1y,
	       const double &B0x, const double &B0y,
	       const double &pre) : density(rho), velocity(vx,vy), pressure(pre),
				    PerturbativeB(B1x,B1y), IntrinsicB(B0x,B0y){ };
  //@}

  //! @name Field access
  //@{
  double& d(void) {return density; }                     //!< Return density variable
  const double& d(void) const {return density; }         //!< Return density (read-only)

  Vector2D& V(void) {return velocity; }			 //!< Return velocity vector
  const Vector2D& V(void) const {return velocity; }	 //!< Return velocity vector (read-only)
  double& vx(void) {return velocity.x; }	         //!< Return x-direction velocity
  const double& vx(void) const {return velocity.x; }	 //!< Return x-direction velocity (read-only)
  double& vy(void) {return velocity.y; }		 //!< Return y-direction velocity
  const double& vy(void) const {return velocity.y; }	 //!< Return y-direction velocity (read-only)

  double& p(void) {return pressure; }                    //!< Return pressure
  const double& p(void) const {return pressure; }        //!< Return pressure (read-only)

  Vector2D& B1(void) {return PerturbativeB; }		 //!< Return perturbative magnetic field
  const Vector2D& B1(void) const {return PerturbativeB;} //!< Return perturbative magnetic field (read-only)
  double& B1x(void) {return PerturbativeB.x; }		 //!< Return x-direction perturbative magnetic field
  const double& B1x(void) const {return PerturbativeB.x;}//!< Return x-direction perturbative magnetic field (read-only)
  double& B1y(void) {return PerturbativeB.y; }		 //!< Return y-direction perturbative magnetic field
  const double& B1y(void) const {return PerturbativeB.y;}//!< Return y-direction perturbative magnetic field (read-only)

  Vector2D& B0(void) {return IntrinsicB; }		 //!< Return intrinsic magnetic field
  const Vector2D& B0(void) const {return IntrinsicB; }	 //!< Return intrinsic magnetic field (read-only)
  double& B0x(void) {return IntrinsicB.x; }		 //!< Return x-direction intrinsic magnetic field
  const double& B0x(void) const {return IntrinsicB.x;}   //!< Return x-direction intrinsic magnetic field (read-only)
  double& B0y(void) {return IntrinsicB.y; }		 //!< Return y-direction intrinsic magnetic field
  const double& B0y(void) const {return IntrinsicB.y;}   //!< Return y-direction intrinsic magnetic field (read-only)
  //@}

  //! Destructor.
  ~MHD2D_pState(void){};

private:

  //! @name Solution state primitive variables:
  //@{
  double          density;       //!< Gas density.
  Vector2D        velocity;      //!< Flow velocity (2D vector).
  double          pressure;      //!< Gas pressure.
  Vector2D        PerturbativeB; //!< Perturbative magnetic field (2D vector).
  Vector2D        IntrinsicB;    //!< Intrinsic magnetic field (2D vector).
  //@}
  
};

#endif  /* _MHD2D_STATE_INCLUDED */

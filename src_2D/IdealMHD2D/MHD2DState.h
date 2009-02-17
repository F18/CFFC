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

#define	NUM_VAR_MHD2D    4

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
 *      B1      -- Return perturbative magnetic field.  
 *      B0      -- Return intrinsic magnetic field.     
 *      p       -- Return pressure.                     
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
  MHD2D_pState(void):density(ONE), velocity(ZERO), PerturbativeB(ZERO), IntrinsicB(ZERO), pressure(ONE){ };

  //! Value Constructor
  explicit MHD2D_pState(const double &Val):density(Val), velocity(Val), PerturbativeB(Val), IntrinsicB(Val), pressure(Val){ };

  //! Copy constructor.
  MHD2D_pState(const MHD2D_pState &W):density(W.density), velocity(W.velocity), PerturbativeB(W.PerturbativeB),
				      IntrinsicB(W.IntrinsicB), pressure(W.pressure){ };
  //@}

  //! Destructor.
  ~MHD2D_pState(void){};

private:
  double          density;       // Gas density.
  Vector2D        velocity;      // Flow velocity (2D vector).
  Vector2D        PerturbativeB; // Perturbative magnetic field (2D vector).
  Vector2D        IntrinsicB;    // Intrinsic magnetic field (2D vector).
  double          pressure;      // Gas pressure.
  
};

#endif  /* _MHD2D_STATE_INCLUDED */

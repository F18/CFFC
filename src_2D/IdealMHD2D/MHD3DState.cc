/*!\file MHD3DState.cc
  @brief Subroutines for 3D Ideal MHD Solution State Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "MHD3DState.h"		// Include 3D ideal MHD solution state header file.


// ======= Member functions/variables of MHD3D_pState class =======

/**************************
 * Assign gas constants.  *
 **************************/
double MHD3D_pState::g = GAMMA_MONATOMIC;
double MHD3D_pState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD3D_pState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);







// ======= Member functions/variables of MHD3D_cState class =======

/**************************
 * Assign gas constants.  *
 **************************/
double MHD3D_cState::g = GAMMA_MONATOMIC;
double MHD3D_cState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD3D_cState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);

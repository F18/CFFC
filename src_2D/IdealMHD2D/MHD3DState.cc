/*!\file MHD3DState.cc
  @brief Subroutines for 3D Ideal MHD Solution State Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "MHD3DState.h"		// Include 3D ideal MHD solution state header file.


/**************************
 * Assign gas constants.  *
 **************************/
double MHD3D_pState::g = GAMMA_MONATOMIC;
double MHD3D_pState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD3D_pState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);


/*!\file MHD2DState.cc
  @brief Subroutines for 2D Ideal MHD Solution State Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "MHD2DState.h"		// Include 2D ideal MHD solution state header file.


/**************************
 * Assign gas constants.  *
 **************************/
double MHD2D_pState::g = GAMMA_MONATOMIC;
double MHD2D_pState::gm1 = GAMMA_MONATOMIC-ONE;
double MHD2D_pState::gm1i = ONE/(GAMMA_MONATOMIC-ONE);


/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Medium2DState.cc                                           **
 **                                                                  **
 ** Description: The radiation state class contains the properties   **
 **              of a gray, absorbing, emitting, anisotropically     **
 **              scattering medium.  This file defines the state     ** 
 **              class.                                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#include "Medium2DState.h"


/********************************************************
 * Static member initialization                         *
 ********************************************************/
int         Medium2D_State :: Nband          = 0;     // number of bands (& quad pts)
int         Medium2D_State :: NUM_VAR_MEDIUM2D = 0;   // total number of variables
SNBCK*      Medium2D_State :: SNBCKdata      = NULL;  // SNBCK data object
double      Medium2D_State :: Absorb_Type    = MEDIUM2D_ABSORB_GRAY; // absorbsion model
Vector2D_Function<Medium2D_State>* Medium2D_State :: Field = NULL;   // functor pointer

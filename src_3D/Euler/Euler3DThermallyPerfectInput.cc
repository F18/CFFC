/****************** Euler3DThermallyPerfectInput.cc *******************
  Constructors for the Euler3D termally perfect gas input class.
***********************************************************************/

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // EULER3D_THERMALLYPERFECT_STATE_INCLUDED

/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
*********************************************************/
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::
Broadcast_Input_Parameters(void) {

   //call the default broadcast input parameters
   Input_Parameters<Euler3D_ThermallyPerfect_pState, 
      Euler3D_ThermallyPerfect_cState>::Broadcast_Input_Parameters();
   
   
#ifdef _MPI_VERSION

   MPI::COMM_WORLD.Bcast(&(Wo.rho), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Wo.v.x), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Wo.v.y), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Wo.v.z), 
                         1, 
                         MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&(Wo.p), 
                         1, 
                         MPI::INT, 0);

#endif

}

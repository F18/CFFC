/****************** Euler3DPolytropicInput.cc *******************
Constructors for the Euler3D termally perfect gas input class.
***********************************************************************/

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED
#include "Euler3DPolytropicState.h"
#endif // EULER3D_POLYTROPIC_STATE_INCLUDED

/********************************************************
* Routine: Broadcast_Input_Parameters                  *
*                                                      *
* Broadcast the input parameters variables to all      *
* processors involved in the calculation from the      *
* primary processor using the MPI broadcast routine.   *
*                                                      *
*********************************************************/
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
Euler3D_Polytropic_cState>::
Broadcast(void) {
	
	//call the default broadcast input parameters
	cout << "\n A1: "; cout.flush();
	cout << "\n CC: "; cout.flush();
	// Input_Parameters<Euler3D_Polytropic_pState, Euler3D_Polytropic_cState>::Broadcast();
	Input_Parameters::Broadcast();
	cout << "\n A2: "; cout.flush();
	
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

/* LES3DFsdHexaBlock.h: Header file defining 3D Hexahedral mesh solution classes
                        for the LES Navier-Stokes equations for a thermally perfect
                        gaseous mixture. */

#ifndef _LES3DFSD_HEXA_BLOCK_INCLUDED
#define _LES3DFSD_HEXA_BLOCK_INCLUDED

// Include required CFFC header files. 

#ifndef _HEXA_INCLUDED
#include "../../HexaBlock/HexaBlock.h"
#endif // HEXA_INCLUDED

#ifndef _LES3DFSD_STATE_INCLUDED
#include "LES3DFsdState.h"
#endif // LES3DFSD_STATE_INCLUDED   

/* Define required specializations. */

template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::allocate_static(void);

template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::deallocate_static(void);

template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
Output_Tecplot(Input_Parameters<LES3DFsd_pState,
      		                LES3DFsd_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3DFsd_pState,
      		                      LES3DFsd_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3DFsd_pState,
      		                      LES3DFsd_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
int Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
ICs(Input_Parameters<LES3DFsd_pState, 
                     LES3DFsd_cState> &IPs);

template<>
int Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
ICs_Specializations(Input_Parameters<LES3DFsd_pState, 
                                     LES3DFsd_cState> &IPs);

template<>
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
BCs(Input_Parameters<LES3DFsd_pState, 
                     LES3DFsd_cState> &IPs);

template<>
double Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
CFL(Input_Parameters<LES3DFsd_pState,
                     LES3DFsd_cState> &IPs);

template<>
int Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3DFsd_pState, 
                                          LES3DFsd_cState> &IPs);

template<>
int Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3DFsd_pState,
                                                     LES3DFsd_cState> &IPs);

template<>
int Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
UnloadReceiveBuffer_Solution(double *buffer,
                             int &buffer_count,
                             const int buffer_size,
                             const int i_min, 
                             const int i_max,
                             const int i_inc,
                             const int j_min, 
                             const int j_max,
                             const int j_inc,
			     const int k_min, 
                             const int k_max,
                             const int k_inc);

#endif /* _LES3D_HEXA_INCLUDED  */

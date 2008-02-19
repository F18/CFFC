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

#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#include "../../TurbulenceModelling/TurbulentVelocityField.h"
#endif // TURBULENT_VELOCITY_FIELD_INCLUDED

#ifndef _TURBULENCE_AVERAGING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceAveraging.h"
#endif // TURBULENCE_AVERAGING_INCLUDED

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


/********************************************************
 *       Burning rate                                   *
 ********************************************************/
template<typename HEXA_BLOCK>
double Turbulent_Burning_Rate(HEXA_BLOCK *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Grid3D_Input_Parameters &IPs){

  double local_vol, Yf_u, rho_u, Ly, Lz, burning_rate(ZERO);
  Yf_u = 0.05518;//Fresh_Fuel_Mass_Fraction;
  rho_u = 1.13;//Fresh_Density;
  Ly = IPs.Box_Width;
  Lz = IPs.Box_Height;
  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     local_vol = Solution_Block[p].Grid.volume(i,j,k);
 	     burning_rate +=  Solution_Block[p].W[i][j][k].Fsd*local_vol*Solution_Block[p].W[i][j][k].rho; 
	   }
	}
      }
    }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = burning_rate*0.3837/(PI*0.0056*0.0056);//(Ly*Lz); //laminar_flame_speed/Ly;//(rho_u*Ly);  //(rho_u*Yf_u*Ly);

  return burning_rate;
}


#endif /* _LES3D_HEXA_INCLUDED  */

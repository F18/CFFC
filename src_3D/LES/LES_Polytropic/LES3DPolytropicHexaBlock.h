/* LES3DPolytropicHexaBlock.h: Header file defining 3D Hexahedral mesh solution classes
                        for the LES Navier-Stokes equations for a polytropic gas. */

#ifndef _LES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED
#define _LES3D_POLYTROPIC_HEXA_BLOCK_INCLUDED

// Include required CFFC header files. 

#ifndef _HEXA_BLOCK_INCLUDED
#include "../../HexaBlock/HexaBlock.h"
#endif // HEXA_BLOCK_INCLUDED

#ifndef _LES3D_POLYTROPIC_STATE_INCLUDED
#include "LES3DPolytropicState.h"
#endif // LES3D_POLYTROPIC_STATE_INCLUDED  

#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#include "../../TurbulenceModelling/TurbulentVelocityField.h"
#endif // TURBULENT_VELOCITY_FIELD_INCLUDED

#ifndef _TURBULENCE_AVERAGING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceAveraging.h"
#endif // TURBULENCE_AVERAGING_INCLUDED

#ifndef _LES_FILTERS_INCLUDED
#include "../Filters/LES_Filters.h"
#endif // TURBULENCE_AVERAGING_INCLUDED

/* Define required specializations. */


template<>
void Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
Output_Tecplot(Input_Parameters<LES3D_Polytropic_pState,
               LES3D_Polytropic_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,  
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File);

template<>
void Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3D_Polytropic_pState,
                     LES3D_Polytropic_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3D_Polytropic_pState,
                     LES3D_Polytropic_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);


template<>
int Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
ICs_Specializations(Input_Parameters<LES3D_Polytropic_pState, 
                                     LES3D_Polytropic_cState> &IPs);

template<>
double Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
CFL(Input_Parameters<LES3D_Polytropic_pState,
                     LES3D_Polytropic_cState> &IPs);

template<>
int Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3D_Polytropic_pState, 
                                          LES3D_Polytropic_cState> &IPs);

template<>
int Hexa_Block<LES3D_Polytropic_pState, LES3D_Polytropic_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3D_Polytropic_pState,
                                                     LES3D_Polytropic_cState> &IPs);



/********************************************************
 *   Return SFS_Kinetic_Energy at specified node        *
 ********************************************************/
template<typename HEXA_BLOCK>
double SFS_Kinetic_Energy_n(HEXA_BLOCK &Soln_Blk,
                     const int &i, 
                     const int &j, 
                     const int &k){
    
    return Trilinear_Interpolation(
                                   Soln_Blk.Grid.Cell[i-1][j][k].Xc, Soln_Blk.W[i-1][j][k].SFS_Kinetic_Energy(Soln_Blk.dWdx[i-1][j][k],
                                                                                                              Soln_Blk.dWdy[i-1][j][k],
                                                                                                              Soln_Blk.dWdz[i-1][j][k],
                                                                                                              Soln_Blk.Grid.Cell[i-1][j][k].V),
                                   Soln_Blk.Grid.Cell[i][j][k].Xc, Soln_Blk.W[i][j][k].SFS_Kinetic_Energy(Soln_Blk.dWdx[i][j][k],
                                                                                                          Soln_Blk.dWdy[i][j][k],
                                                                                                          Soln_Blk.dWdz[i][j][k],
                                                                                                          Soln_Blk.Grid.Cell[i][j][k].V),
                                   Soln_Blk.Grid.Cell[i][j-1][k].Xc, Soln_Blk.W[i][j-1][k].SFS_Kinetic_Energy(Soln_Blk.dWdx[i][j-1][k],
                                                                                                              Soln_Blk.dWdy[i][j-1][k],
                                                                                                              Soln_Blk.dWdz[i][j-1][k],
                                                                                                              Soln_Blk.Grid.Cell[i][j-1][k].V),
                                   Soln_Blk.Grid.Cell[i-1][j-1][k].Xc, Soln_Blk.W[i-1][j-1][k].SFS_Kinetic_Energy(Soln_Blk.dWdx[i-1][j-1][k],
                                                                                                                  Soln_Blk.dWdy[i-1][j-1][k],
                                                                                                                  Soln_Blk.dWdz[i-1][j-1][k],
                                                                                                                  Soln_Blk.Grid.Cell[i-1][j-1][k].V),
                                   Soln_Blk.Grid.Cell[i-1][j][k-1].Xc, Soln_Blk.W[i-1][j][k-1].SFS_Kinetic_Energy(Soln_Blk.dWdx[i-1][j][k-1],
                                                                                                                  Soln_Blk.dWdy[i-1][j][k-1],
                                                                                                                  Soln_Blk.dWdz[i-1][j][k-1],
                                                                                                                  Soln_Blk.Grid.Cell[i-1][j][k-1].V),
                                   Soln_Blk.Grid.Cell[i][j][k-1].Xc, Soln_Blk.W[i][j][k-1].SFS_Kinetic_Energy(Soln_Blk.dWdx[i][j][k-1],
                                                                                                              Soln_Blk.dWdy[i][j][k-1],
                                                                                                              Soln_Blk.dWdz[i][j][k-1],
                                                                                                              Soln_Blk.Grid.Cell[i][j][k-1].V),
                                   Soln_Blk.Grid.Cell[i][j-1][k-1].Xc, Soln_Blk.W[i][j-1][k-1].SFS_Kinetic_Energy(Soln_Blk.dWdx[i][j-1][k-1],
                                                                                                                  Soln_Blk.dWdy[i][j-1][k-1],
                                                                                                                  Soln_Blk.dWdz[i][j-1][k-1],
                                                                                                                  Soln_Blk.Grid.Cell[i][j-1][k-1].V),
                                   Soln_Blk.Grid.Cell[i-1][j-1][k-1].Xc, Soln_Blk.W[i-1][j-1][k-1].SFS_Kinetic_Energy(Soln_Blk.dWdx[i-1][j-1][k-1],
                                                                                                                      Soln_Blk.dWdy[i-1][j-1][k-1],
                                                                                                                      Soln_Blk.dWdz[i-1][j-1][k-1],
                                                                                                                      Soln_Blk.Grid.Cell[i-1][j-1][k-1].V),
                                   Soln_Blk.Grid.Node[i][j][k].X);
    
}



#endif /* _LES3D_HEXA_INCLUDED  */

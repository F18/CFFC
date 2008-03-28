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

class HexaBlock_Node : public Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> {

  public:

   double Reaction_Rate_Fsd_n(const int i,
                              const int j,
                              const int k);

   double M_x_n(const int i,
                const int j,
                const int k);

   double M_y_n(const int i,
                const int j,
                const int k);

   double M_z_n(const int i,
                const int j,
                const int k);

   double Resolved_Strain_n(const int i,
 		            const int j,
			    const int k);

   double Resolved_Propagation_Curvature_n(const int i,
                                           const int j,
                                           const int k);

   double SFS_Strain_n(const int i,
                       const int j,
                       const int k,
                       const int Flow_Type);

   double SFS_Curvature_n(const int i,
                          const int j,
                          const int k);

   double Resolved_Curvature_n(const int i,
                               const int j,
                               const int k);

   double Resolved_Propagation_n(const int i,
                                 const int j,
                                 const int k);

   double Resolved_Convection_Progvar_n(const int i,
                                        const int j,
                                        const int k);   

   double Resolved_Convection_Fsd_n(const int i,
                                    const int j,
                                    const int k);

   double NGT_Progvar_n(const int i,
                        const int j,
                        const int k);

   double NGT_Fsd_n(const int i,
                    const int j,
                    const int k);

   double SFS_Diffusion_Progvar_n(const int i,
                                  const int j,
                                  const int k,
                                  const int Flow_Type);

   double SFS_Diffusion_Fsd_n(const int i,
                              const int j,
                              const int k,
                              const int Flow_Type);

   double Heat_Release_Strain_n(const int i,
                                const int j,
                                const int k);

   double Net_Rate_Change_Progvar_n(const int i,
                                    const int j,
                                    const int k,
                                    const int Flow_Type);

   double Net_Rate_Change_Fsd_n(const int i,
                                const int j,
                                const int k,
                                const int Flow_Type);

   //constructors

   HexaBlock_Node(void){ Hexa_Block(); }
   HexaBlock_Node(Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> &Soln_Blk){ Copy_static(Soln_Blk); }

   //destructors

   ~HexaBlock_Node(void){ deallocate(); deallocate_static(); }
 
};

/********************************************************
 *        Return Q_criterion at specified node          *
 ********************************************************/
template<typename HEXA_BLOCK>
double Q_criterion_n(HEXA_BLOCK &Soln_Blk,
		     const int &i, 
		     const int &j, 
		     const int &k){

  return Trilinear_Interpolation(
    Soln_Blk.Grid.Cell[i-1][j][k].Xc, Soln_Blk.W[i-1][j][k].Q_criterion(Soln_Blk.dWdx[i-1][j][k],
									Soln_Blk.dWdy[i-1][j][k],
									Soln_Blk.dWdz[i-1][j][k]),
    Soln_Blk.Grid.Cell[i][j][k].Xc, Soln_Blk.W[i][j][k].Q_criterion(Soln_Blk.dWdx[i][j][k],
								    Soln_Blk.dWdy[i][j][k],
								    Soln_Blk.dWdz[i][j][k]),
    Soln_Blk.Grid.Cell[i][j-1][k].Xc, Soln_Blk.W[i][j-1][k].Q_criterion(Soln_Blk.dWdx[i][j-1][k],
									Soln_Blk.dWdy[i][j-1][k],
									Soln_Blk.dWdz[i][j-1][k]),
    Soln_Blk.Grid.Cell[i-1][j-1][k].Xc, Soln_Blk.W[i-1][j-1][k].Q_criterion(Soln_Blk.dWdx[i-1][j-1][k],
									    Soln_Blk.dWdy[i-1][j-1][k],
									    Soln_Blk.dWdz[i-1][j-1][k]),
    Soln_Blk.Grid.Cell[i-1][j][k-1].Xc, Soln_Blk.W[i-1][j][k-1].Q_criterion(Soln_Blk.dWdx[i-1][j][k-1],
									    Soln_Blk.dWdy[i-1][j][k-1],
									    Soln_Blk.dWdz[i-1][j][k-1]),
    Soln_Blk.Grid.Cell[i][j][k-1].Xc, Soln_Blk.W[i][j][k-1].Q_criterion(Soln_Blk.dWdx[i][j][k-1],
									Soln_Blk.dWdy[i][j][k-1],
									Soln_Blk.dWdz[i][j][k-1]),
    Soln_Blk.Grid.Cell[i][j-1][k-1].Xc, Soln_Blk.W[i][j-1][k-1].Q_criterion(Soln_Blk.dWdx[i][j-1][k-1],
									    Soln_Blk.dWdy[i][j-1][k-1],
									    Soln_Blk.dWdz[i][j-1][k-1]),
    Soln_Blk.Grid.Cell[i-1][j-1][k-1].Xc, Soln_Blk.W[i-1][j-1][k-1].Q_criterion(Soln_Blk.dWdx[i-1][j-1][k-1],
										Soln_Blk.dWdy[i-1][j-1][k-1],
										Soln_Blk.dWdz[i-1][j-1][k-1]),
    Soln_Blk.Grid.Node[i][j][k].X);

}

/********************************************************
 *  Return Reaction_Rate_Fsd at specified node          *
 ********************************************************/
inline double HexaBlock_Node::Reaction_Rate_Fsd_n(const int i, 
						  const int j, 
						  const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Reaction_Rate_Fsd(dWdx[i-1][j][k],
								dWdy[i-1][j][k],
								dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Reaction_Rate_Fsd(dWdx[i][j][k],
								dWdy[i][j][k],
								dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Reaction_Rate_Fsd(dWdx[i][j-1][k],
								dWdy[i][j-1][k],
								dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Reaction_Rate_Fsd(dWdx[i-1][j-1][k],
								dWdy[i-1][j-1][k],
								dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Reaction_Rate_Fsd(dWdx[i-1][j][k-1],
								dWdy[i-1][j][k-1],
								dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Reaction_Rate_Fsd(dWdx[i][j][k-1],
								dWdy[i][j][k-1],
								dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Reaction_Rate_Fsd(dWdx[i][j-1][k-1],
									dWdy[i][j-1][k-1],
									dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Reaction_Rate_Fsd(dWdx[i-1][j-1][k-1],
									dWdy[i-1][j-1][k-1],
									dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *                Return M_x at specified node          *
 ********************************************************/
inline double HexaBlock_Node::M_x_n(const int i,
				    const int j,
				    const int k){
  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].M_x(dWdx[i-1][j][k],
					      dWdy[i-1][j][k],
					      dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].M_x(dWdx[i][j][k],
					  dWdy[i][j][k],
					  dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].M_x(dWdx[i][j-1][k],
					      dWdy[i][j-1][k],
					      dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].M_x(dWdx[i-1][j-1][k],
						  dWdy[i-1][j-1][k],
						  dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].M_x(dWdx[i-1][j][k-1],
						  dWdy[i-1][j][k-1],
						  dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].M_x(dWdx[i][j][k-1],
					      dWdy[i][j][k-1],
					      dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].M_x(dWdx[i][j-1][k-1],
						  dWdy[i][j-1][k-1],
						  dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].M_x(dWdx[i-1][j-1][k-1],
						      dWdy[i-1][j-1][k-1],
						      dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *                Return M_y at specified node          *
 ********************************************************/
inline double HexaBlock_Node::M_y_n(const int i,
				    const int j,
				    const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].M_y(dWdx[i-1][j][k],
					      dWdy[i-1][j][k],
					      dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].M_y(dWdx[i][j][k],
					  dWdy[i][j][k],
					  dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].M_y(dWdx[i][j-1][k],
					      dWdy[i][j-1][k],
					      dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].M_y(dWdx[i-1][j-1][k],
						  dWdy[i-1][j-1][k],
						  dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].M_y(dWdx[i-1][j][k-1],
						  dWdy[i-1][j][k-1],
						  dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].M_y(dWdx[i][j][k-1],
					      dWdy[i][j][k-1],
					      dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].M_y(dWdx[i][j-1][k-1],
						  dWdy[i][j-1][k-1],
						  dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].M_y(dWdx[i-1][j-1][k-1],
						      dWdy[i-1][j-1][k-1],
						      dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *                Return M_z at specified node          *
 ********************************************************/
inline double HexaBlock_Node::M_z_n(const int i,
				    const int j,
				    const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].M_z(dWdx[i-1][j][k],
					      dWdy[i-1][j][k],
					      dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].M_z(dWdx[i][j][k],
					  dWdy[i][j][k],
					  dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].M_z(dWdx[i][j-1][k],
					      dWdy[i][j-1][k],
					      dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].M_z(dWdx[i-1][j-1][k],
						  dWdy[i-1][j-1][k],
						  dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].M_z(dWdx[i-1][j][k-1],
						  dWdy[i-1][j][k-1],
						  dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].M_z(dWdx[i][j][k-1],
					      dWdy[i][j][k-1],
					      dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].M_z(dWdx[i][j-1][k-1],
						  dWdy[i][j-1][k-1],
						  dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].M_z(dWdx[i-1][j-1][k-1],
						      dWdy[i-1][j-1][k-1],
						      dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *    Return Resolved_Strain at specified node          *
 ********************************************************/
inline double HexaBlock_Node::Resolved_Strain_n(const int i,
						const int j,
						const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Strain(dWdx[i-1][j][k],
							  dWdy[i-1][j][k],
							  dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Strain(dWdx[i][j][k],
						      dWdy[i][j][k],
						      dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Strain(dWdx[i][j-1][k],
							  dWdy[i][j-1][k],
							  dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Strain(dWdx[i-1][j-1][k],
							      dWdy[i-1][j-1][k],
							      dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Strain(dWdx[i-1][j][k-1],
							      dWdy[i-1][j][k-1],
							      dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Strain(dWdx[i][j][k-1],
							  dWdy[i][j][k-1],
							  dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Strain(dWdx[i][j-1][k-1],
							      dWdy[i][j-1][k-1],
							      dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Strain(dWdx[i-1][j-1][k-1],
								  dWdy[i-1][j-1][k-1],
								  dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************************
 * Return Resolved_Propagation_Curvature at specified node          *
 ********************************************************************/
inline double HexaBlock_Node::Resolved_Propagation_Curvature_n(const int i,
							       const int j,
							       const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Propagation_Curvature(dWdx[i-1][j][k],
									 dWdy[i-1][j][k],
									 dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Propagation_Curvature(dWdx[i][j][k],
								     dWdy[i][j][k],
								     dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Propagation_Curvature(dWdx[i][j-1][k],
									 dWdy[i][j-1][k],
									 dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Propagation_Curvature(dWdx[i-1][j-1][k],
									     dWdy[i-1][j-1][k],
									     dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Propagation_Curvature(dWdx[i-1][j][k-1],
									     dWdy[i-1][j][k-1],
									     dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Propagation_Curvature(dWdx[i][j][k-1],
									 dWdy[i][j][k-1],
									 dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Propagation_Curvature(dWdx[i][j-1][k-1],
									     dWdy[i][j-1][k-1],
									     dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Propagation_Curvature(dWdx[i-1][j-1][k-1],
										 dWdy[i-1][j-1][k-1],
										 dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *        Return SFS_Strain at specified node           *
 ********************************************************/
inline double HexaBlock_Node::SFS_Strain_n(const int i,
					   const int j,
					   const int k,
					   const int Flow_Type){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].SFS_Strain(dWdx[i-1][j][k],
						     dWdy[i-1][j][k],
						     dWdz[i-1][j][k],
						     Flow_Type,
						     Grid.volume(i-1,j,k)),
    Grid.Cell[i][j][k].Xc, W[i][j][k].SFS_Strain(dWdx[i][j][k],
						 dWdy[i][j][k],
						 dWdz[i][j][k],
						 Flow_Type,
						 Grid.volume(i,j,k)),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].SFS_Strain(dWdx[i][j-1][k],
						     dWdy[i][j-1][k],
						     dWdz[i][j-1][k],
						     Flow_Type,
						     Grid.volume(i,j-1,k)),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].SFS_Strain(dWdx[i-1][j-1][k],
							 dWdy[i-1][j-1][k],
							 dWdz[i-1][j-1][k],
							 Flow_Type,
							 Grid.volume(i-1,j-1,k)),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].SFS_Strain(dWdx[i-1][j][k-1],
							 dWdy[i-1][j][k-1],
							 dWdz[i-1][j][k-1],
							 Flow_Type,
							 Grid.volume(i-1,j,k-1)),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].SFS_Strain(dWdx[i][j][k-1],
						     dWdy[i][j][k-1],
						     dWdz[i][j][k-1],
						     Flow_Type,
						     Grid.volume(i,j,k-1)),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].SFS_Strain(dWdx[i][j-1][k-1],
							 dWdy[i][j-1][k-1],
							 dWdz[i][j-1][k-1],
							 Flow_Type,
							 Grid.volume(i,j-1,k-1)),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].SFS_Strain(dWdx[i-1][j-1][k-1],
							     dWdy[i-1][j-1][k-1],
							     dWdz[i-1][j-1][k-1],
							     Flow_Type,
							     Grid.volume(i-1,j-1,k-1)),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *      Return SFS_Curvature at specified node          *
 ********************************************************/
inline double HexaBlock_Node::SFS_Curvature_n(const int i,
					      const int j,
					      const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].SFS_Curvature(dWdx[i-1][j][k],
							dWdy[i-1][j][k],
							dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].SFS_Curvature(dWdx[i][j][k],
						    dWdy[i][j][k],
						    dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].SFS_Curvature(dWdx[i][j-1][k],
							dWdy[i][j-1][k],
							dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].SFS_Curvature(dWdx[i-1][j-1][k],
							    dWdy[i-1][j-1][k],
							    dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].SFS_Curvature(dWdx[i-1][j][k-1],
							    dWdy[i-1][j][k-1],
							    dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].SFS_Curvature(dWdx[i][j][k-1],
							dWdy[i][j][k-1],
							dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].SFS_Curvature(dWdx[i][j-1][k-1],
							    dWdy[i][j-1][k-1],
							    dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].SFS_Curvature(dWdx[i-1][j-1][k-1],
								dWdy[i-1][j-1][k-1],
								dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/***************************************************************
 *        Return Resolved_Curvature at specified node          *
 ***************************************************************/
inline double HexaBlock_Node::Resolved_Curvature_n(const int i,
						   const int j,
						   const int k){
  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Curvature(dWdx[i-1][j][k],
							     dWdy[i-1][j][k],
							     dWdz[i-1][j][k],
							     _d2Wdx2[i-1][j][k],
							     _d2Wdy2[i-1][j][k],
							     _d2Wdz2[i-1][j][k],
							     _d2Wdxdy[i-1][j][k],
							     _d2Wdxdz[i-1][j][k],
							     _d2Wdydz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Curvature(dWdx[i][j][k],
							 dWdy[i][j][k],
							 dWdz[i][j][k],
							 _d2Wdx2[i][j][k],
							 _d2Wdy2[i][j][k],
							 _d2Wdz2[i][j][k],
							 _d2Wdxdy[i][j][k],
							 _d2Wdxdz[i][j][k],
							 _d2Wdydz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Curvature(dWdx[i][j-1][k],
							     dWdy[i][j-1][k],
							     dWdz[i][j-1][k],
							     _d2Wdx2[i][j-1][k],
							     _d2Wdy2[i][j-1][k],
							     _d2Wdz2[i][j-1][k],
							     _d2Wdxdy[i][j-1][k],
							     _d2Wdxdz[i][j-1][k],
							     _d2Wdydz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Curvature(dWdx[i-1][j-1][k],
								 dWdy[i-1][j-1][k],
								 dWdz[i-1][j-1][k],
								 _d2Wdx2[i-1][j-1][k],
								 _d2Wdy2[i-1][j-1][k],
								 _d2Wdz2[i-1][j-1][k],
								 _d2Wdxdy[i-1][j-1][k],
								 _d2Wdxdz[i-1][j-1][k],
								 _d2Wdydz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Curvature(dWdx[i-1][j][k-1],
								 dWdy[i-1][j][k-1],
								 dWdz[i-1][j][k-1],
								 _d2Wdx2[i-1][j][k-1],
								 _d2Wdy2[i-1][j][k-1],
								 _d2Wdz2[i-1][j][k-1],
								 _d2Wdxdy[i-1][j][k-1],
								 _d2Wdxdz[i-1][j][k-1],
								 _d2Wdydz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Curvature(dWdx[i][j][k-1],
							     dWdy[i][j][k-1],
							     dWdz[i][j][k-1],
							     _d2Wdx2[i][j][k-1],
							     _d2Wdy2[i][j][k-1],
							     _d2Wdz2[i][j][k-1],
							     _d2Wdxdy[i][j][k-1],
							     _d2Wdxdz[i][j][k-1],
							     _d2Wdydz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Curvature(dWdx[i][j-1][k-1],
								 dWdy[i][j-1][k-1],
								 dWdz[i][j-1][k-1],
								 _d2Wdx2[i][j-1][k-1],
								 _d2Wdy2[i][j-1][k-1],
								 _d2Wdz2[i][j-1][k-1],
								 _d2Wdxdy[i][j-1][k-1],
								 _d2Wdxdz[i][j-1][k-1],
								 _d2Wdydz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Curvature(dWdx[i-1][j-1][k-1],
								     dWdy[i-1][j-1][k-1],
								     dWdz[i-1][j-1][k-1],
								     _d2Wdx2[i-1][j-1][k-1],
								     _d2Wdy2[i-1][j-1][k-1],
								     _d2Wdz2[i-1][j-1][k-1],
								     _d2Wdxdy[i-1][j-1][k-1],
								     _d2Wdxdz[i-1][j-1][k-1],
								     _d2Wdydz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/***************************************************************
 *      Return Resolved_Propagation at specified node          *
 ***************************************************************/
inline double HexaBlock_Node::Resolved_Propagation_n(const int i,
						     const int j,
						     const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Propagation(dWdx[i-1][j][k],
							       dWdy[i-1][j][k],
							       dWdz[i-1][j][k],
							       _d2Wdx2[i-1][j][k],
							       _d2Wdy2[i-1][j][k],
							       _d2Wdz2[i-1][j][k],
							       _d2Wdxdy[i-1][j][k],
							       _d2Wdxdz[i-1][j][k],
							       _d2Wdydz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Propagation(dWdx[i][j][k],
							   dWdy[i][j][k],
							   dWdz[i][j][k],
							   _d2Wdx2[i][j][k],
							   _d2Wdy2[i][j][k],
							   _d2Wdz2[i][j][k],
							   _d2Wdxdy[i][j][k],
							   _d2Wdxdz[i][j][k],
							   _d2Wdydz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Propagation(dWdx[i][j-1][k],
							       dWdy[i][j-1][k],
							       dWdz[i][j-1][k],
							       _d2Wdx2[i][j-1][k],
							       _d2Wdy2[i][j-1][k],
							       _d2Wdz2[i][j-1][k],
							       _d2Wdxdy[i][j-1][k],
							       _d2Wdxdz[i][j-1][k],
							       _d2Wdydz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Propagation(dWdx[i-1][j-1][k],
								   dWdy[i-1][j-1][k],
								   dWdz[i-1][j-1][k],
								   _d2Wdx2[i-1][j-1][k],
								   _d2Wdy2[i-1][j-1][k],
								   _d2Wdz2[i-1][j-1][k],
								   _d2Wdxdy[i-1][j-1][k],
								   _d2Wdxdz[i-1][j-1][k],
								   _d2Wdydz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Propagation(dWdx[i-1][j][k-1],
								   dWdy[i-1][j][k-1],
								   dWdz[i-1][j][k-1],
								   _d2Wdx2[i-1][j][k-1],
								   _d2Wdy2[i-1][j][k-1],
								   _d2Wdz2[i-1][j][k-1],
								   _d2Wdxdy[i-1][j][k-1],
								   _d2Wdxdz[i-1][j][k-1],
								   _d2Wdydz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Propagation(dWdx[i][j][k-1],
							       dWdy[i][j][k-1],
							       dWdz[i][j][k-1],
							       _d2Wdx2[i][j][k-1],
							       _d2Wdy2[i][j][k-1],
							       _d2Wdz2[i][j][k-1],
							       _d2Wdxdy[i][j][k-1],
							       _d2Wdxdz[i][j][k-1],
							       _d2Wdydz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Propagation(dWdx[i][j-1][k-1],
								   dWdy[i][j-1][k-1],
								   dWdz[i][j-1][k-1],
								   _d2Wdx2[i][j-1][k-1],
								   _d2Wdy2[i][j-1][k-1],
								   _d2Wdz2[i][j-1][k-1],
								   _d2Wdxdy[i][j-1][k-1],
								   _d2Wdxdz[i][j-1][k-1],
								   _d2Wdydz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Propagation(dWdx[i-1][j-1][k-1],
								       dWdy[i-1][j-1][k-1],
								       dWdz[i-1][j-1][k-1],
								       _d2Wdx2[i-1][j-1][k-1],
								       _d2Wdy2[i-1][j-1][k-1],
								       _d2Wdz2[i-1][j-1][k-1],
								       _d2Wdxdy[i-1][j-1][k-1],
								       _d2Wdxdz[i-1][j-1][k-1],
								       _d2Wdydz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/************************************************************************
 *        Return Resolved_Convection_Progvar at specified node          *
 ************************************************************************/
inline double HexaBlock_Node::Resolved_Convection_Progvar_n(const int i,
							    const int j,
							    const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Convection_Progvar(dWdx[i-1][j][k],
								      dWdy[i-1][j][k],
								      dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Convection_Progvar(dWdx[i][j][k],
								  dWdy[i][j][k],
								  dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Convection_Progvar(dWdx[i][j-1][k],
								      dWdy[i][j-1][k],
								      dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Convection_Progvar(dWdx[i-1][j-1][k],
									  dWdy[i-1][j-1][k],
									  dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Convection_Progvar(dWdx[i-1][j][k-1],
									  dWdy[i-1][j][k-1],
									  dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Convection_Progvar(dWdx[i][j][k-1],
								      dWdy[i][j][k-1],
								      dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Convection_Progvar(dWdx[i][j-1][k-1],
									  dWdy[i][j-1][k-1],
									  dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Convection_Progvar(dWdx[i-1][j-1][k-1],
									      dWdy[i-1][j-1][k-1],
									      dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}
/********************************************************************
 *        Return Resolved_Convection_Fsd at specified node          *
 ********************************************************************/
inline double HexaBlock_Node::Resolved_Convection_Fsd_n(const int i,
							const int j,
							const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Resolved_Convection_Fsd(dWdx[i-1][j][k],
								  dWdy[i-1][j][k],
								  dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Resolved_Convection_Fsd(dWdx[i][j][k],
							      dWdy[i][j][k],
							      dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Resolved_Convection_Fsd(dWdx[i][j-1][k],
								  dWdy[i][j-1][k],
								  dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Resolved_Convection_Fsd(dWdx[i-1][j-1][k],
								      dWdy[i-1][j-1][k],
								      dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Resolved_Convection_Fsd(dWdx[i-1][j][k-1],
								      dWdy[i-1][j][k-1],
								      dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Resolved_Convection_Fsd(dWdx[i][j][k-1],
								  dWdy[i][j][k-1],
								  dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Resolved_Convection_Fsd(dWdx[i][j-1][k-1],
								      dWdy[i][j-1][k-1],
								      dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Resolved_Convection_Fsd(dWdx[i-1][j-1][k-1],
									  dWdy[i-1][j-1][k-1],
									  dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}
/*******************************************************************
 *                   Return NGT_Progvar at specified node          *
 *******************************************************************/
inline double HexaBlock_Node::NGT_Progvar_n(const int i,
					    const int j,
					    const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].NGT_Progvar(dWdx[i-1][j][k],
						      dWdy[i-1][j][k],
						      dWdz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].NGT_Progvar(dWdx[i][j][k],
						  dWdy[i][j][k],
						  dWdz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].NGT_Progvar(dWdx[i][j-1][k],
						      dWdy[i][j-1][k],
						      dWdz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].NGT_Progvar(dWdx[i-1][j-1][k],
							  dWdy[i-1][j-1][k],
							  dWdz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].NGT_Progvar(dWdx[i-1][j][k-1],
							  dWdy[i-1][j][k-1],
							  dWdz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].NGT_Progvar(dWdx[i][j][k-1],
						      dWdy[i][j][k-1],
						      dWdz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].NGT_Progvar(dWdx[i][j-1][k-1],
							  dWdy[i][j-1][k-1],
							  dWdz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].NGT_Progvar(dWdx[i-1][j-1][k-1],
							      dWdy[i-1][j-1][k-1],
							      dWdz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/***************************************************************
 *                   Return NGT_Fsd at specified node          *
 ***************************************************************/
inline double HexaBlock_Node::NGT_Fsd_n(const int i,
					const int j,
					const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].NGT_Fsd(dWdx[i-1][j][k],
						  dWdy[i-1][j][k],
						  dWdz[i-1][j][k],
						  _d2Wdx2[i-1][j][k],
						  _d2Wdy2[i-1][j][k],
						  _d2Wdz2[i-1][j][k],
						  _d2Wdxdy[i-1][j][k],
						  _d2Wdxdz[i-1][j][k],
						  _d2Wdydz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].NGT_Fsd(dWdx[i][j][k],
					      dWdy[i][j][k],
					      dWdz[i][j][k],
					      _d2Wdx2[i][j][k],
					      _d2Wdy2[i][j][k],
					      _d2Wdz2[i][j][k],
					      _d2Wdxdy[i][j][k],
					      _d2Wdxdz[i][j][k],
					      _d2Wdydz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].NGT_Fsd(dWdx[i][j-1][k],
						  dWdy[i][j-1][k],
						  dWdz[i][j-1][k],
						  _d2Wdx2[i][j-1][k],
						  _d2Wdy2[i][j-1][k],
						  _d2Wdz2[i][j-1][k],
						  _d2Wdxdy[i][j-1][k],
						  _d2Wdxdz[i][j-1][k],
						  _d2Wdydz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].NGT_Fsd(dWdx[i-1][j-1][k],
						      dWdy[i-1][j-1][k],
						      dWdz[i-1][j-1][k],
						      _d2Wdx2[i-1][j-1][k],
						      _d2Wdy2[i-1][j-1][k],
						      _d2Wdz2[i-1][j-1][k],
						      _d2Wdxdy[i-1][j-1][k],
						      _d2Wdxdz[i-1][j-1][k],
						      _d2Wdydz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].NGT_Fsd(dWdx[i-1][j][k-1],
						      dWdy[i-1][j][k-1],
						      dWdz[i-1][j][k-1],
						      _d2Wdx2[i-1][j][k-1],
						      _d2Wdy2[i-1][j][k-1],
						      _d2Wdz2[i-1][j][k-1],
						      _d2Wdxdy[i-1][j][k-1],
						      _d2Wdxdz[i-1][j][k-1],
						      _d2Wdydz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].NGT_Fsd(dWdx[i][j][k-1],
						  dWdy[i][j][k-1],
						  dWdz[i][j][k-1],
						  _d2Wdx2[i][j][k-1],
						  _d2Wdy2[i][j][k-1],
						  _d2Wdz2[i][j][k-1],
						  _d2Wdxdy[i][j][k-1],
						  _d2Wdxdz[i][j][k-1],
						  _d2Wdydz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].NGT_Fsd(dWdx[i][j-1][k-1],
						      dWdy[i][j-1][k-1],
						      dWdz[i][j-1][k-1],
						      _d2Wdx2[i][j-1][k-1],
						      _d2Wdy2[i][j-1][k-1],
						      _d2Wdz2[i][j-1][k-1],
						      _d2Wdxdy[i][j-1][k-1],
						      _d2Wdxdz[i][j-1][k-1],
						      _d2Wdydz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].NGT_Fsd(dWdx[i-1][j-1][k-1],
							  dWdy[i-1][j-1][k-1],
							  dWdz[i-1][j-1][k-1],
							  _d2Wdx2[i-1][j-1][k-1],
							  _d2Wdy2[i-1][j-1][k-1],
							  _d2Wdz2[i-1][j-1][k-1],
							  _d2Wdxdy[i-1][j-1][k-1],
							  _d2Wdxdz[i-1][j-1][k-1],
							  _d2Wdydz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/*****************************************************************************
 *                   Return SFS_Diffusion_Progvar at specified node          *
 *****************************************************************************/
inline double HexaBlock_Node::SFS_Diffusion_Progvar_n(const int i,
						      const int j,
						      const int k,
						      const int Flow_Typs){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].SFS_Diffusion_Progvar(dWdx[i-1][j][k],
								dWdy[i-1][j][k],
								dWdz[i-1][j][k],
								_d2Wdx2[i-1][j][k],
								_d2Wdy2[i-1][j][k],
								_d2Wdz2[i-1][j][k],
								_d2Wdxdy[i-1][j][k],
								_d2Wdxdz[i-1][j][k],
								_d2Wdydz[i-1][j][k],
								Flow_Type,
								Grid.volume(i-1,j,k)),
    Grid.Cell[i][j][k].Xc, W[i][j][k].SFS_Diffusion_Progvar(dWdx[i][j][k],
							    dWdy[i][j][k],
							    dWdz[i][j][k],
							    _d2Wdx2[i][j][k],
							    _d2Wdy2[i][j][k],
							    _d2Wdz2[i][j][k],
							    _d2Wdxdy[i][j][k],
							    _d2Wdxdz[i][j][k],
							    _d2Wdydz[i][j][k],
							    Flow_Type,
							    Grid.volume(i,j,k)),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].SFS_Diffusion_Progvar(dWdx[i][j-1][k],
								dWdy[i][j-1][k],
								dWdz[i][j-1][k],
								_d2Wdx2[i][j-1][k],
								_d2Wdy2[i][j-1][k],
								_d2Wdz2[i][j-1][k],
								_d2Wdxdy[i][j-1][k],
								_d2Wdxdz[i][j-1][k],
								_d2Wdydz[i][j-1][k],
								Flow_Type,
								Grid.volume(i,j-1,k)),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].SFS_Diffusion_Progvar(dWdx[i-1][j-1][k],
								    dWdy[i-1][j-1][k],
								    dWdz[i-1][j-1][k],
								    _d2Wdx2[i-1][j-1][k],
								    _d2Wdy2[i-1][j-1][k],
								    _d2Wdz2[i-1][j-1][k],
								    _d2Wdxdy[i-1][j-1][k],
								    _d2Wdxdz[i-1][j-1][k],
								    _d2Wdydz[i-1][j-1][k],
								    Flow_Type,
								    Grid.volume(i-1,j-1,k)),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].SFS_Diffusion_Progvar(dWdx[i-1][j][k-1],
								    dWdy[i-1][j][k-1],
								    dWdz[i-1][j][k-1],
								    _d2Wdx2[i-1][j][k-1],
								    _d2Wdy2[i-1][j][k-1],
								    _d2Wdz2[i-1][j][k-1],
								    _d2Wdxdy[i-1][j][k-1],
								    _d2Wdxdz[i-1][j][k-1],
								    _d2Wdydz[i-1][j][k-1],
								    Flow_Type,
								    Grid.volume(i-1,j,k-1)),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].SFS_Diffusion_Progvar(dWdx[i][j][k-1],
								dWdy[i][j][k-1],
								dWdz[i][j][k-1],
								_d2Wdx2[i][j][k-1],
								_d2Wdy2[i][j][k-1],
								_d2Wdz2[i][j][k-1],
								_d2Wdxdy[i][j][k-1],
								_d2Wdxdz[i][j][k-1],
								_d2Wdydz[i][j][k-1],
								Flow_Type,
								Grid.volume(i,j,k-1)),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].SFS_Diffusion_Progvar(dWdx[i][j-1][k-1],
								    dWdy[i][j-1][k-1],
								    dWdz[i][j-1][k-1],
								    _d2Wdx2[i][j-1][k-1],
								    _d2Wdy2[i][j-1][k-1],
								    _d2Wdz2[i][j-1][k-1],
								    _d2Wdxdy[i][j-1][k-1],
								    _d2Wdxdz[i][j-1][k-1],
								    _d2Wdydz[i][j-1][k-1],
								    Flow_Type,
								    Grid.volume(i,j-1,k-1)),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].SFS_Diffusion_Progvar(dWdx[i-1][j-1][k-1],
									dWdy[i-1][j-1][k-1],
									dWdz[i-1][j-1][k-1],
									_d2Wdx2[i-1][j-1][k-1],
									_d2Wdy2[i-1][j-1][k-1],
									_d2Wdz2[i-1][j-1][k-1],
									_d2Wdxdy[i-1][j-1][k-1],
									_d2Wdxdz[i-1][j-1][k-1],
									_d2Wdydz[i-1][j-1][k-1],
									Flow_Type,
									Grid.volume(i-1,j-1,k-1)),
    Grid.Node[i][j][k].X);

}

/*****************************************************************************
 *                         Return SFS_Diffusion_Fsd at specified node        *
 *****************************************************************************/
inline double HexaBlock_Node::SFS_Diffusion_Fsd_n(const int i,
						  const int j,
						  const int k,
						  const int Flow_Typs){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].SFS_Diffusion_Fsd(dWdx[i-1][j][k],
							    dWdy[i-1][j][k],
							    dWdz[i-1][j][k],
							    _d2Wdx2[i-1][j][k],
							    _d2Wdy2[i-1][j][k],
							    _d2Wdz2[i-1][j][k],
							    _d2Wdxdy[i-1][j][k],
							    _d2Wdxdz[i-1][j][k],
							    _d2Wdydz[i-1][j][k],
							    Flow_Type,
							    Grid.volume(i-1,j,k)),
    Grid.Cell[i][j][k].Xc, W[i][j][k].SFS_Diffusion_Fsd(dWdx[i][j][k],
							dWdy[i][j][k],
							dWdz[i][j][k],
							_d2Wdx2[i][j][k],
							_d2Wdy2[i][j][k],
							_d2Wdz2[i][j][k],
							_d2Wdxdy[i][j][k],
							_d2Wdxdz[i][j][k],
							_d2Wdydz[i][j][k],
							Flow_Type,
							Grid.volume(i,j,k)),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].SFS_Diffusion_Fsd(dWdx[i][j-1][k],
							    dWdy[i][j-1][k],
							    dWdz[i][j-1][k],
							    _d2Wdx2[i][j-1][k],
							    _d2Wdy2[i][j-1][k],
							    _d2Wdz2[i][j-1][k],
							    _d2Wdxdy[i][j-1][k],
							    _d2Wdxdz[i][j-1][k],
							    _d2Wdydz[i][j-1][k],
							    Flow_Type,
							    Grid.volume(i,j-1,k)),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].SFS_Diffusion_Fsd(dWdx[i-1][j-1][k],
								dWdy[i-1][j-1][k],
								dWdz[i-1][j-1][k],
								_d2Wdx2[i-1][j-1][k],
								_d2Wdy2[i-1][j-1][k],
								_d2Wdz2[i-1][j-1][k],
								_d2Wdxdy[i-1][j-1][k],
								_d2Wdxdz[i-1][j-1][k],
								_d2Wdydz[i-1][j-1][k],
								Flow_Type,
								Grid.volume(i-1,j-1,k)),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].SFS_Diffusion_Fsd(dWdx[i-1][j][k-1],
								dWdy[i-1][j][k-1],
								dWdz[i-1][j][k-1],
								_d2Wdx2[i-1][j][k-1],
								_d2Wdy2[i-1][j][k-1],
								_d2Wdz2[i-1][j][k-1],
								_d2Wdxdy[i-1][j][k-1],
								_d2Wdxdz[i-1][j][k-1],
								_d2Wdydz[i-1][j][k-1],
								Flow_Type,
								Grid.volume(i-1,j,k-1)),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].SFS_Diffusion_Fsd(dWdx[i][j][k-1],
							    dWdy[i][j][k-1],
							    dWdz[i][j][k-1],
							    _d2Wdx2[i][j][k-1],
							    _d2Wdy2[i][j][k-1],
							    _d2Wdz2[i][j][k-1],
							    _d2Wdxdy[i][j][k-1],
							    _d2Wdxdz[i][j][k-1],
							    _d2Wdydz[i][j][k-1],
							    Flow_Type,
							    Grid.volume(i,j,k-1)),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].SFS_Diffusion_Fsd(dWdx[i][j-1][k-1],
								dWdy[i][j-1][k-1],
								dWdz[i][j-1][k-1],
								_d2Wdx2[i][j-1][k-1],
								_d2Wdy2[i][j-1][k-1],
								_d2Wdz2[i][j-1][k-1],
								_d2Wdxdy[i][j-1][k-1],
								_d2Wdxdz[i][j-1][k-1],
								_d2Wdydz[i][j-1][k-1],
								Flow_Type,
								Grid.volume(i,j-1,k-1)),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].SFS_Diffusion_Fsd(dWdx[i-1][j-1][k-1],
								    dWdy[i-1][j-1][k-1],
								    dWdz[i-1][j-1][k-1],
								    _d2Wdx2[i-1][j-1][k-1],
								    _d2Wdy2[i-1][j-1][k-1],
								    _d2Wdz2[i-1][j-1][k-1],
								    _d2Wdxdy[i-1][j-1][k-1],
								    _d2Wdxdz[i-1][j-1][k-1],
								    _d2Wdydz[i-1][j-1][k-1],
								    Flow_Type,
								    Grid.volume(i-1,j-1,k-1)),
    Grid.Node[i][j][k].X);

}

/***************************************************************************
 *                   Return Heat_Release_Strain at specified node          *
 ***************************************************************************/
inline double HexaBlock_Node::Heat_Release_Strain_n(const int i,
						    const int j,
						    const int k){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Heat_Release_Strain(dWdx[i-1][j][k],
							      dWdy[i-1][j][k],
							      dWdz[i-1][j][k],
							      _d2Wdx2[i-1][j][k],
							      _d2Wdy2[i-1][j][k],
							      _d2Wdz2[i-1][j][k],
							      _d2Wdxdy[i-1][j][k],
							      _d2Wdxdz[i-1][j][k],
							      _d2Wdydz[i-1][j][k]),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Heat_Release_Strain(dWdx[i][j][k],
							  dWdy[i][j][k],
							  dWdz[i][j][k],
							  _d2Wdx2[i][j][k],
							  _d2Wdy2[i][j][k],
							  _d2Wdz2[i][j][k],
							  _d2Wdxdy[i][j][k],
							  _d2Wdxdz[i][j][k],
							  _d2Wdydz[i][j][k]),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Heat_Release_Strain(dWdx[i][j-1][k],
							      dWdy[i][j-1][k],
							      dWdz[i][j-1][k],
							      _d2Wdx2[i][j-1][k],
							      _d2Wdy2[i][j-1][k],
							      _d2Wdz2[i][j-1][k],
							      _d2Wdxdy[i][j-1][k],
							      _d2Wdxdz[i][j-1][k],
							      _d2Wdydz[i][j-1][k]),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Heat_Release_Strain(dWdx[i-1][j-1][k],
								  dWdy[i-1][j-1][k],
								  dWdz[i-1][j-1][k],
								  _d2Wdx2[i-1][j-1][k],
								  _d2Wdy2[i-1][j-1][k],
								  _d2Wdz2[i-1][j-1][k],
								  _d2Wdxdy[i-1][j-1][k],
								  _d2Wdxdz[i-1][j-1][k],
								  _d2Wdydz[i-1][j-1][k]),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Heat_Release_Strain(dWdx[i-1][j][k-1],
								  dWdy[i-1][j][k-1],
								  dWdz[i-1][j][k-1],
								  _d2Wdx2[i-1][j][k-1],
								  _d2Wdy2[i-1][j][k-1],
								  _d2Wdz2[i-1][j][k-1],
								  _d2Wdxdy[i-1][j][k-1],
								  _d2Wdxdz[i-1][j][k-1],
								  _d2Wdydz[i-1][j][k-1]),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Heat_Release_Strain(dWdx[i][j][k-1],
							      dWdy[i][j][k-1],
							      dWdz[i][j][k-1],
							      _d2Wdx2[i][j][k-1],
							      _d2Wdy2[i][j][k-1],
							      _d2Wdz2[i][j][k-1],
							      _d2Wdxdy[i][j][k-1],
							      _d2Wdxdz[i][j][k-1],
							      _d2Wdydz[i][j][k-1]),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Heat_Release_Strain(dWdx[i][j-1][k-1],
								  dWdy[i][j-1][k-1],
								  dWdz[i][j-1][k-1],
								  _d2Wdx2[i][j-1][k-1],
								  _d2Wdy2[i][j-1][k-1],
								  _d2Wdz2[i][j-1][k-1],
								  _d2Wdxdy[i][j-1][k-1],
								  _d2Wdxdz[i][j-1][k-1],
								  _d2Wdydz[i][j-1][k-1]),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Heat_Release_Strain(dWdx[i-1][j-1][k-1],
								      dWdy[i-1][j-1][k-1],
								      dWdz[i-1][j-1][k-1],
								      _d2Wdx2[i-1][j-1][k-1],
								      _d2Wdy2[i-1][j-1][k-1],
								      _d2Wdz2[i-1][j-1][k-1],
								      _d2Wdxdy[i-1][j-1][k-1],
								      _d2Wdxdz[i-1][j-1][k-1],
								      _d2Wdydz[i-1][j-1][k-1]),
    Grid.Node[i][j][k].X);

}

/*****************************************************************************
 *                   Return Net_Rate_Change_Progvar at specified node          *
 *****************************************************************************/
inline double HexaBlock_Node::Net_Rate_Change_Progvar_n(const int i,
							const int j,
							const int k,
							const int Flow_Typs){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Net_Rate_Change_Progvar(dWdx[i-1][j][k],
								  dWdy[i-1][j][k],
								  dWdz[i-1][j][k],
								  _d2Wdx2[i-1][j][k],
								  _d2Wdy2[i-1][j][k],
								  _d2Wdz2[i-1][j][k],
								  _d2Wdxdy[i-1][j][k],
								  _d2Wdxdz[i-1][j][k],
								  _d2Wdydz[i-1][j][k],
								  Flow_Type,
								  Grid.volume(i-1,j,k)),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Net_Rate_Change_Progvar(dWdx[i][j][k],
							      dWdy[i][j][k],
							      dWdz[i][j][k],
							      _d2Wdx2[i][j][k],
							      _d2Wdy2[i][j][k],
							      _d2Wdz2[i][j][k],
							      _d2Wdxdy[i][j][k],
							      _d2Wdxdz[i][j][k],
							      _d2Wdydz[i][j][k],
							      Flow_Type,
							      Grid.volume(i,j,k)),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Net_Rate_Change_Progvar(dWdx[i][j-1][k],
								  dWdy[i][j-1][k],
								  dWdz[i][j-1][k],
								  _d2Wdx2[i][j-1][k],
								  _d2Wdy2[i][j-1][k],
								  _d2Wdz2[i][j-1][k],
								  _d2Wdxdy[i][j-1][k],
								  _d2Wdxdz[i][j-1][k],
								  _d2Wdydz[i][j-1][k],
								  Flow_Type,
								  Grid.volume(i,j-1,k)),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Net_Rate_Change_Progvar(dWdx[i-1][j-1][k],
								      dWdy[i-1][j-1][k],
								      dWdz[i-1][j-1][k],
								      _d2Wdx2[i-1][j-1][k],
								      _d2Wdy2[i-1][j-1][k],
								      _d2Wdz2[i-1][j-1][k],
								      _d2Wdxdy[i-1][j-1][k],
								      _d2Wdxdz[i-1][j-1][k],
								      _d2Wdydz[i-1][j-1][k],
								      Flow_Type,
								      Grid.volume(i-1,j-1,k)),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Net_Rate_Change_Progvar(dWdx[i-1][j][k-1],
								      dWdy[i-1][j][k-1],
								      dWdz[i-1][j][k-1],
								      _d2Wdx2[i-1][j][k-1],
								      _d2Wdy2[i-1][j][k-1],
								      _d2Wdz2[i-1][j][k-1],
								      _d2Wdxdy[i-1][j][k-1],
								      _d2Wdxdz[i-1][j][k-1],
								      _d2Wdydz[i-1][j][k-1],
								      Flow_Type,
								      Grid.volume(i-1,j,k-1)),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Net_Rate_Change_Progvar(dWdx[i][j][k-1],
								  dWdy[i][j][k-1],
								  dWdz[i][j][k-1],
								  _d2Wdx2[i][j][k-1],
								  _d2Wdy2[i][j][k-1],
								  _d2Wdz2[i][j][k-1],
								  _d2Wdxdy[i][j][k-1],
								  _d2Wdxdz[i][j][k-1],
								  _d2Wdydz[i][j][k-1],
								  Flow_Type,
								  Grid.volume(i,j,k-1)),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Net_Rate_Change_Progvar(dWdx[i][j-1][k-1],
								      dWdy[i][j-1][k-1],
								      dWdz[i][j-1][k-1],
								      _d2Wdx2[i][j-1][k-1],
								      _d2Wdy2[i][j-1][k-1],
								      _d2Wdz2[i][j-1][k-1],
								      _d2Wdxdy[i][j-1][k-1],
								      _d2Wdxdz[i][j-1][k-1],
								      _d2Wdydz[i][j-1][k-1],
								      Flow_Type,
								      Grid.volume(i,j-1,k-1)),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Net_Rate_Change_Progvar(dWdx[i-1][j-1][k-1],
									  dWdy[i-1][j-1][k-1],
									  dWdz[i-1][j-1][k-1],
									  _d2Wdx2[i-1][j-1][k-1],
									  _d2Wdy2[i-1][j-1][k-1],
									  _d2Wdz2[i-1][j-1][k-1],
									  _d2Wdxdy[i-1][j-1][k-1],
									  _d2Wdxdz[i-1][j-1][k-1],
									  _d2Wdydz[i-1][j-1][k-1],
									  Flow_Type,
									  Grid.volume(i-1,j-1,k-1)),
    Grid.Node[i][j][k].X);

}

/*****************************************************************************
 *                     Return Net_Rate_Change_Fsd at specified node          *
 *****************************************************************************/
inline double HexaBlock_Node::Net_Rate_Change_Fsd_n(const int i,
						    const int j,
						    const int k,
						    const int Flow_Typs){

  return Trilinear_Interpolation(
    Grid.Cell[i-1][j][k].Xc, W[i-1][j][k].Net_Rate_Change_Fsd(dWdx[i-1][j][k],
							      dWdy[i-1][j][k],
							      dWdz[i-1][j][k],
							      _d2Wdx2[i-1][j][k],
							      _d2Wdy2[i-1][j][k],
							      _d2Wdz2[i-1][j][k],
							      _d2Wdxdy[i-1][j][k],
							      _d2Wdxdz[i-1][j][k],
							      _d2Wdydz[i-1][j][k],
							      Flow_Type,
							      Grid.volume(i-1,j,k)),
    Grid.Cell[i][j][k].Xc, W[i][j][k].Net_Rate_Change_Fsd(dWdx[i][j][k],
							  dWdy[i][j][k],
							  dWdz[i][j][k],
							  _d2Wdx2[i][j][k],
							  _d2Wdy2[i][j][k],
							  _d2Wdz2[i][j][k],
							  _d2Wdxdy[i][j][k],
							  _d2Wdxdz[i][j][k],
							  _d2Wdydz[i][j][k],
							  Flow_Type,
							  Grid.volume(i,j,k)),
    Grid.Cell[i][j-1][k].Xc, W[i][j-1][k].Net_Rate_Change_Fsd(dWdx[i][j-1][k],
							      dWdy[i][j-1][k],
							      dWdz[i][j-1][k],
							      _d2Wdx2[i][j-1][k],
							      _d2Wdy2[i][j-1][k],
							      _d2Wdz2[i][j-1][k],
							      _d2Wdxdy[i][j-1][k],
							      _d2Wdxdz[i][j-1][k],
							      _d2Wdydz[i][j-1][k],
							      Flow_Type,
							      Grid.volume(i,j-1,k)),
    Grid.Cell[i-1][j-1][k].Xc, W[i-1][j-1][k].Net_Rate_Change_Fsd(dWdx[i-1][j-1][k],
								  dWdy[i-1][j-1][k],
								  dWdz[i-1][j-1][k],
								  _d2Wdx2[i-1][j-1][k],
								  _d2Wdy2[i-1][j-1][k],
								  _d2Wdz2[i-1][j-1][k],
								  _d2Wdxdy[i-1][j-1][k],
								  _d2Wdxdz[i-1][j-1][k],
								  _d2Wdydz[i-1][j-1][k],
								  Flow_Type,
								  Grid.volume(i-1,j-1,k)),
    Grid.Cell[i-1][j][k-1].Xc, W[i-1][j][k-1].Net_Rate_Change_Fsd(dWdx[i-1][j][k-1],
								  dWdy[i-1][j][k-1],
								  dWdz[i-1][j][k-1],
								  _d2Wdx2[i-1][j][k-1],
								  _d2Wdy2[i-1][j][k-1],
								  _d2Wdz2[i-1][j][k-1],
								  _d2Wdxdy[i-1][j][k-1],
								  _d2Wdxdz[i-1][j][k-1],
								  _d2Wdydz[i-1][j][k-1],
								  Flow_Type,
								  Grid.volume(i-1,j,k-1)),
    Grid.Cell[i][j][k-1].Xc, W[i][j][k-1].Net_Rate_Change_Fsd(dWdx[i][j][k-1],
							      dWdy[i][j][k-1],
							      dWdz[i][j][k-1],
							      _d2Wdx2[i][j][k-1],
							      _d2Wdy2[i][j][k-1],
							      _d2Wdz2[i][j][k-1],
							      _d2Wdxdy[i][j][k-1],
							      _d2Wdxdz[i][j][k-1],
							      _d2Wdydz[i][j][k-1],
							      Flow_Type,
							      Grid.volume(i,j,k-1)),
    Grid.Cell[i][j-1][k-1].Xc, W[i][j-1][k-1].Net_Rate_Change_Fsd(dWdx[i][j-1][k-1],
								  dWdy[i][j-1][k-1],
								  dWdz[i][j-1][k-1],
								  _d2Wdx2[i][j-1][k-1],
								  _d2Wdy2[i][j-1][k-1],
								  _d2Wdz2[i][j-1][k-1],
								  _d2Wdxdy[i][j-1][k-1],
								  _d2Wdxdz[i][j-1][k-1],
								  _d2Wdydz[i][j-1][k-1],
								  Flow_Type,
								  Grid.volume(i,j-1,k-1)),
    Grid.Cell[i-1][j-1][k-1].Xc, W[i-1][j-1][k-1].Net_Rate_Change_Fsd(dWdx[i-1][j-1][k-1],
								      dWdy[i-1][j-1][k-1],
								      dWdz[i-1][j-1][k-1],
								      _d2Wdx2[i-1][j-1][k-1],
								      _d2Wdy2[i-1][j-1][k-1],
								      _d2Wdz2[i-1][j-1][k-1],
								      _d2Wdxdy[i-1][j-1][k-1],
								      _d2Wdxdz[i-1][j-1][k-1],
								      _d2Wdydz[i-1][j-1][k-1],
								      Flow_Type,
								      Grid.volume(i-1,j-1,k-1)),
    Grid.Node[i][j][k].X);

}

/********************************************************
 *      propagation_dir_area                            *
 ********************************************************/
template<typename HEXA_BLOCK>
double propagation_dir_area(HEXA_BLOCK &Solution_Block,
			    const int &i,
			    const int &j,
			    const int &k) {

  double area, temp_dot;

  Vector3D N_iso_Yfuel(Solution_Block.dWdx[i][j][k].spec[0].c,
		       Solution_Block.dWdy[i][j][k].spec[0].c,
		       Solution_Block.dWdz[i][j][k].spec[0].c);

  // East face  
  double dot_prod = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceE(i,j,k));
  area = Solution_Block.Grid.AfaceE(i,j,k);

  // West face
  temp_dot = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceW(i,j,k));
  if ( temp_dot > dot_prod ) { 
    dot_prod = temp_dot;
    area = Solution_Block.Grid.AfaceW(i,j,k);
  }

  // North face
  temp_dot = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceN(i,j,k));
  if ( temp_dot > dot_prod ) { 
    dot_prod = temp_dot;
    area = Solution_Block.Grid.AfaceN(i,j,k);
  }

  // South face
  temp_dot = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceS(i,j,k));
  if ( temp_dot > dot_prod ) { 
    dot_prod = temp_dot;
    area = Solution_Block.Grid.AfaceS(i,j,k);
  }

  // Top face
  temp_dot = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceTop(i,j,k));
  if ( temp_dot > dot_prod ) { 
    dot_prod = temp_dot;
    area = Solution_Block.Grid.AfaceTop(i,j,k);
  }

  // Bottom
  temp_dot = N_iso_Yfuel.dot(Solution_Block.Grid.nfaceBot(i,j,k));
  if ( temp_dot > dot_prod ) { 
    dot_prod = temp_dot;
    area = Solution_Block.Grid.AfaceBot(i,j,k);
  }

  return area;
}

/********************************************************
 *       Burning rate                                   *
 ********************************************************/
template<typename HEXA_BLOCK>
double Turbulent_Burning_Rate(HEXA_BLOCK *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Grid3D_Input_Parameters &IPs){

  double local_vol, Yf_u, rho_u, burning_rate(ZERO), iso_surface_area(ZERO);
  double flame_height = 0.004;
  Yf_u = 0.05518;//Fresh_Fuel_Mass_Fraction;
  rho_u = 1.13;//Fresh_Density;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     local_vol = Solution_Block[p].Grid.volume(i,j,k);
 	     burning_rate +=  Solution_Block[p].W[i][j][k].Fsd*Solution_Block[p].W[i][j][k].rho*local_vol; 
	     if (Solution_Block[p].W[i][j][k].C <= 0.5 && 
                 Solution_Block[p].Grid.Cell[i][j][k].Xc.z > flame_height) {
	         flame_height = Solution_Block[p].Grid.Cell[i][j][k].Xc.z;
/* 	       iso_surface_area += propagation_dir_area(Solution_Block[p], i, j, k); */
	     }
	   }
	}
      }
    }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = burning_rate*0.403/(PI*0.0056*(0.0056+2.0*flame_height));

/*   iso_surface_area = CFFC_Summation_MPI(iso_surface_area); */
  
/*   double ref_area, Lx, Ly, Lz; */

/*   if ( IPs.i_Grid == GRID_PERIODIC_BOX_WITH_INFLOW  || */
/*        IPs.i_Grid == GRID_PERIODIC_BOX) { */
/*     Ly = IPs.Box_Height; */
/*     Lz = IPs.Box_Length; */
/*     ref_area = Ly*Lz; */
/*   } else if ( IPs.i_Grid == GRID_BUNSEN_BOX ) { */
/*     Lx = IPs.Box_Width; */
/*     ref_area = (0.025 + 2.0*0.02)*Lx; */
/*   } else if ( IPs.i_Grid == GRID_BUNSEN_BURNER ) { */
/*     ref_area = PI*0.0056*0.0056*0.005; */
/*   } */

/*   if ( iso_surface_area > ref_area ) { */
/*     burning_rate = -burning_rate/(rho_u*Yf_u*iso_surface_area); */
/*   } else { */
/*     burning_rate = -burning_rate/(rho_u*Yf_u*ref_area); */
/*   } */

  return burning_rate;
}


#endif /* _LES3D_HEXA_INCLUDED  */

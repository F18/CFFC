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
void Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>::
Linear_Reconstruction_LeastSquares(const int i,
				   const int j,
				   const int k,
				   const int Limiter);

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

   HexaBlock_Node(void){ Hexa_Block<LES3DFsd_pState, LES3DFsd_cState>(); }
   HexaBlock_Node(Hexa_Block<LES3DFsd_pState, LES3DFsd_cState> &Soln_Blk){ Copy_static(Soln_Blk); }

   //destructors

   ~HexaBlock_Node(void){ deallocate(); deallocate_static(); }
 
};

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


#endif /* _LES3D_HEXA_INCLUDED  */

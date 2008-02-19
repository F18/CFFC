/* LES3DThickenedFlameHexaBlock.h: Header file defining 3D Hexahedral mesh solution classes
                                   for the LES Navier-Stokes equations for a thermally perfect
                                   gaseous mixture. */

#ifndef _LES3DTF_HEXA_BLOCK_INCLUDED
#define _LES3DTF_HEXA_BLOCK_INCLUDED

// Include required CFFC header files. 

#ifndef _HEXA_INCLUDED
#include "../../HexaBlock/HexaBlock.h"
#endif // HEXA_INCLUDED

#ifndef _LES3DTF_STATE_INCLUDED
#include "LES3DThickenedFlameState.h"
#endif // LES3DTF_STATE_INCLUDED   

#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#include "../../TurbulenceModelling/TurbulentVelocityField.h"
#endif // TURBULENT_VELOCITY_FIELD_INCLUDED

#ifndef _TURBULENCE_AVERAGING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceAveraging.h"
#endif // TURBULENCE_AVERAGING_INCLUDED

/* Define required specializations. */

template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Output_Tecplot(Input_Parameters<LES3DTF_pState,
      		                LES3DTF_cState> &IPs,
               const int Number_of_Time_Steps,
	       const double &Time,  
               const int Block_Number,
	       const int Output_Title,
	       ostream &Out_File);
				
template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3DTF_pState,
      		                      LES3DTF_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3DTF_pState,
      		                      LES3DTF_cState> &IPs, 
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File);

template<>
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
ICs(Input_Parameters<LES3DTF_pState, 
                     LES3DTF_cState> &IPs);

template<>
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
ICs_Specializations(Input_Parameters<LES3DTF_pState, 
                                     LES3DTF_cState> &IPs);

template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
BCs(Input_Parameters<LES3DTF_pState, 
                     LES3DTF_cState> &IPs);

template<>
double Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
CFL(Input_Parameters<LES3DTF_pState,
                     LES3DTF_cState> &IPs);

template<>
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3DTF_pState, 
                                          LES3DTF_cState> &IPs);

template<>
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3DTF_pState,
                                                     LES3DTF_cState> &IPs);

template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Linear_Reconstruction_LeastSquares(const int i,
				   const int j,
				   const int k,
				   const int Limiter);

/******************************************************************
 * Hexa_Block::NumVar -- Return the number of solution variables. *
 ******************************************************************/
template<>
int Hexa_Block<LES3DTF_pState,LES3DTF_cState>::NumVar(void){
  return (W[0][0][0].num_vars+2);
}



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
 *       Magnitude of the Laplacian of vorticity        *
 *               using finite difference                * 
 ********************************************************/
template<typename HEXA_BLOCK>
double Laplacian_of_Vorticity_FD(HEXA_BLOCK &Solution_Block,
				 const int &i, 
				 const int &j, 
				 const int &k) {
   
  Vector3D vorticity_plus, vorticity_minus, vorticity_local, Laplacian_vorticity;
  double DX,DY,DZ;
  double d2vorx_dxx, d2vory_dxx, d2vorz_dxx;
  double d2vorx_dyy, d2vory_dyy, d2vorz_dyy;
  double d2vorx_dzz, d2vory_dzz, d2vorz_dzz;

  
  vorticity_local = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k],
							Solution_Block.dWdy[i][j][k],
							Solution_Block.dWdz[i][j][k]);
  // d2_dxx
  vorticity_plus = Solution_Block.W[i+1][j][k].vorticity(Solution_Block.dWdx[i+1][j][k],
							 Solution_Block.dWdy[i+1][j][k],
							 Solution_Block.dWdz[i+1][j][k]);
  vorticity_minus = Solution_Block.W[i-1][j][k].vorticity(Solution_Block.dWdx[i-1][j][k],
							  Solution_Block.dWdy[i-1][j][k],
							  Solution_Block.dWdz[i-1][j][k]);  
  DX = 0.5*(Solution_Block.Grid.Cell[i+1][j][k].Xc.x - Solution_Block.Grid.Cell[i-1][j][k].Xc.x);
  d2vorx_dxx = (vorticity_plus.x - 2.0*vorticity_local.x + vorticity_minus.x)/(DX*DX);
  d2vory_dxx = (vorticity_plus.y - 2.0*vorticity_local.y + vorticity_minus.y)/(DX*DX);
  d2vorz_dxx = (vorticity_plus.z - 2.0*vorticity_local.z + vorticity_minus.z)/(DX*DX);

  // d2_dyy
  vorticity_plus = Solution_Block.W[i][j+1][k].vorticity(Solution_Block.dWdx[i][j+1][k],
							 Solution_Block.dWdy[i][j+1][k],
							 Solution_Block.dWdz[i][j+1][k]);
  vorticity_minus = Solution_Block.W[i][j-1][k].vorticity(Solution_Block.dWdx[i][j-1][k],
							  Solution_Block.dWdy[i][j-1][k],
							  Solution_Block.dWdz[i][j-1][k]);  
  DY = 0.5*(Solution_Block.Grid.Cell[i][j+1][k].Xc.y - Solution_Block.Grid.Cell[i][j-1][k].Xc.y);
  d2vorx_dyy = (vorticity_plus.x - 2.0*vorticity_local.x + vorticity_minus.x)/(DY*DY);
  d2vory_dyy = (vorticity_plus.y - 2.0*vorticity_local.y + vorticity_minus.y)/(DY*DY);
  d2vorz_dyy = (vorticity_plus.z - 2.0*vorticity_local.z + vorticity_minus.z)/(DY*DY);

  // d2_dzz
  vorticity_plus = Solution_Block.W[i][j][k+1].vorticity(Solution_Block.dWdx[i][j][k+1],
							 Solution_Block.dWdy[i][j][k+1],
							 Solution_Block.dWdz[i][j][k+1]);
  vorticity_minus = Solution_Block.W[i][j][k-1].vorticity(Solution_Block.dWdx[i][j][k-1],
							  Solution_Block.dWdy[i][j][k-1],
							  Solution_Block.dWdz[i][j][k-1]);  
  DZ = 0.5*(Solution_Block.Grid.Cell[i][j][k+1].Xc.z - Solution_Block.Grid.Cell[i][j][k-1].Xc.z);
  d2vorx_dzz = (vorticity_plus.x - 2.0*vorticity_local.x + vorticity_minus.x)/(DZ*DZ);
  d2vory_dzz = (vorticity_plus.y - 2.0*vorticity_local.y + vorticity_minus.y)/(DZ*DZ);
  d2vorz_dzz = (vorticity_plus.z - 2.0*vorticity_local.z + vorticity_minus.z)/(DZ*DZ);
  

  // components of Laplacian
  Laplacian_vorticity.x = d2vorx_dxx + d2vorx_dyy + d2vorx_dzz;
  Laplacian_vorticity.y = d2vory_dxx + d2vory_dyy + d2vory_dzz;
  Laplacian_vorticity.z = d2vorz_dxx + d2vorz_dyy + d2vorz_dzz;
 

  return Laplacian_vorticity.abs();
  
}


/********************************************************
 *       Magnitude of the Laplacian of vorticity        *
 ********************************************************/
template<typename HEXA_BLOCK>
double Laplacian_of_Vorticity(HEXA_BLOCK &Solution_Block,
			      const int &i,
			      const int &j,
			      const int &k); 

/********************************************************
 *    Gradient of vorticity using least squares         *
 ********************************************************/
template<typename HEXA_BLOCK>
void Gradient_of_Vorticity(const HEXA_BLOCK &Solution_Block,
			   const int &i,
			   const int &j,
			   const int &k,
			   Vector3D &dVdx,
			   Vector3D &dVdy,
			   Vector3D &dVdz);


/********************************************************
 *       Gradient at a face of a cell                   *
 ********************************************************/
template<typename T>
void Face_Gradient(const T &Wc,
		   const T &Wc_Neighbour,
		   const T &dWdx,
		   const T &dWdy,
		   const T &dWdz,
		   const T &dWdx_Neighbour,
		   const T &dWdy_Neighbour,
		   const T &dWdz_Neighbour,
		   const Vector3D &norm,
		   const Vector3D &ts, 
		   const double &deltad, 
		   const double &Volume, 
		   const double &Volume_Neighbour,
		   T &dWdx_face,
		   T &dWdy_face,
		   T &dWdz_face) {

  // weighted factor based on volume
  double alpha = Volume/(Volume + Volume_Neighbour);
   
  Vector3D dWdx_Weighted, dWdy_Weighted, dWdz_Weighted, Grad_middle_term;
  
  dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neighbour;
  dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neighbour;
  dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neighbour;
   
  // Evaluate a weighted term for solution gradients  
  Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y + dWdz_Weighted*ts.z;
      
  // Evaluate gradients of primitive variables on the face
  dWdx_face = (Wc_Neighbour - Wc)/deltad *norm.x/dot(norm, ts) + 
              (dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
  dWdy_face = (Wc_Neighbour - Wc)/deltad *norm.y/dot(norm, ts) + 
              (dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
  dWdz_face = (Wc_Neighbour - Wc)/deltad *norm.z/dot(norm, ts) + 
              (dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));

}


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
	     burning_rate += Solution_Block[p].W[i][j][k].Sw(Solution_Block[p].W[i][j][k].React.reactset_flag).rhospec[0].c*local_vol;
	   }
	}
      }
    }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  burning_rate = -burning_rate/(rho_u*Yf_u*Ly*Lz);

  return burning_rate;

}


#endif /* _LES3D_HEXA_INCLUDED  */

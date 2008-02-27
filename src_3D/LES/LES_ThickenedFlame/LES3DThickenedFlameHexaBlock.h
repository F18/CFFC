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
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
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


template<class HEXA_BLOCK>
double Laplacian_of_Vorticity(HEXA_BLOCK &Solution_Block,
			      const int &i, 
			      const int &j, 
			      const int &k) {
   
  Vector3D vorticity_plus, vorticity_minus, vorticity_local, Laplacian_vorticity;
  double DX,DY,DZ;

/*   double dvorx_dx_plus, dvory_dx_plus, dvorz_dx_plus; */
/*   double dvorx_dx_minus, dvory_dx_minus, dvorz_dx_minus; */

/*   double dvorx_dy_plus, dvory_dy_plus, dvorz_dy_plus; */
/*   double dvorx_dy_minus, dvory_dy_minus, dvorz_dy_minus; */

/*   double dvorx_dz_plus, dvory_dz_plus, dvorz_dz_plus; */
/*   double dvorx_dz_minus, dvory_dz_minus, dvorz_dz_minus; */

  double d2vorx_dxx, d2vory_dxx, d2vorz_dxx;
  double d2vorx_dyy, d2vory_dyy, d2vorz_dyy;
  double d2vorx_dzz, d2vory_dzz, d2vorz_dzz;

/*   // plus  x */
/*   vorticity_plus = Solution_Block.W[i+1][j][k].vorticity(Solution_Block.dWdx[i+1][j][k], */
/* 							 Solution_Block.dWdy[i+1][j][k], */
/* 							 Solution_Block.dWdz[i+1][j][k]); */
/*   vorticity_minus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], */
/* 							Solution_Block.dWdy[i][j][k], */
/* 							Solution_Block.dWdz[i][j][k]); */
/*   DX = Solution_Block.Grid.Cell[i+1][j][k].Xc.x - Solution_Block.Grid.Cell[i][j][k].Xc.x; */
/*   dvorx_dx_plus = (vorticity_plus.x - vorticity_minus.x)/DX; */
/*   dvory_dx_plus = (vorticity_plus.y - vorticity_minus.y)/DX; */
/*   dvorz_dx_plus = (vorticity_plus.z - vorticity_minus.z)/DX; */


/*   // plus y */
/*   vorticity_plus = Solution_Block.W[i][j+1][k].vorticity(Solution_Block.dWdx[i][j+1][k], */
/* 							 Solution_Block.dWdy[i][j+1][k], */
/* 							 Solution_Block.dWdz[i][j+1][k]); */
/*   /\* vorticity_minus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], *\/ */
/* /\* 							Solution_Block.dWdy[i][j][k], *\/ */
/* /\* 							Solution_Block.dWdz[i][j][k]); *\/ */
/*   DY = Solution_Block.Grid.Cell[i][j+1][k].Xc.y - Solution_Block.Grid.Cell[i][j][k].Xc.y; */
/*   dvorx_dy_plus = (vorticity_plus.x - vorticity_minus.x)/DY; */
/*   dvory_dy_plus = (vorticity_plus.y - vorticity_minus.y)/DY; */
/*   dvorz_dy_plus = (vorticity_plus.z - vorticity_minus.z)/DY; */


/*   // plus z */
/*   vorticity_plus = Solution_Block.W[i][j][k+1].vorticity(Solution_Block.dWdx[i][j][k+1], */
/* 							 Solution_Block.dWdy[i][j][k+1], */
/* 							 Solution_Block.dWdz[i][j][k+1]); */
/*   /\* vorticity_minus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], *\/ */
/* /\* 							Solution_Block.dWdy[i][j][k], *\/ */
/* /\* 							Solution_Block.dWdz[i][j][k]); *\/ */
/*   DZ = Solution_Block.Grid.Cell[i][j][k+1].Xc.z - Solution_Block.Grid.Cell[i][j][k].Xc.z; */
/*   dvorx_dz_plus = (vorticity_plus.x - vorticity_minus.x)/DZ; */
/*   dvory_dz_plus = (vorticity_plus.y - vorticity_minus.y)/DZ; */
/*   dvorz_dz_plus = (vorticity_plus.z - vorticity_minus.z)/DZ; */


 
/*   // minus x */
/*   vorticity_plus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], */
/* 						       Solution_Block.dWdy[i][j][k], */
/* 						       Solution_Block.dWdz[i][j][k]); */
/*   vorticity_minus = Solution_Block.W[i-1][j][k].vorticity(Solution_Block.dWdx[i-1][j][k], */
/* 							  Solution_Block.dWdy[i-1][j][k], */
/* 							  Solution_Block.dWdz[i-1][j][k]); */
/*   DX = Solution_Block.Grid.Cell[i][j][k].Xc.x - Solution_Block.Grid.Cell[i-1][j][k].Xc.x; */
/*   dvorx_dx_minus = (vorticity_plus.x - vorticity_minus.x)/DX; */
/*   dvory_dx_minus = (vorticity_plus.y - vorticity_minus.y)/DX; */
/*   dvorz_dx_minus = (vorticity_plus.z - vorticity_minus.z)/DX; */


/*   // minus y */
/*   /\* vorticity_plus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], *\/ */
/* /\* 						       Solution_Block.dWdy[i][j][k], *\/ */
/* /\* 						       Solution_Block.dWdz[i][j][k]); *\/ */
/*   vorticity_minus = Solution_Block.W[i][j-1][k].vorticity(Solution_Block.dWdx[i][j-1][k], */
/* 							  Solution_Block.dWdy[i][j-1][k], */
/* 							  Solution_Block.dWdz[i][j-1][k]); */
/*   DY = Solution_Block.Grid.Cell[i][j][k].Xc.y - Solution_Block.Grid.Cell[i][j-1][k].Xc.y; */
/*   dvorx_dy_minus = (vorticity_plus.x - vorticity_minus.x)/DY; */
/*   dvory_dy_minus = (vorticity_plus.y - vorticity_minus.y)/DY; */
/*   dvorz_dy_minus = (vorticity_plus.z - vorticity_minus.z)/DY; */


/*   // minus z */
/*   /\* vorticity_plus = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k], *\/ */
/* /\* 						       Solution_Block.dWdy[i][j][k], *\/ */
/* /\* 						       Solution_Block.dWdz[i][j][k]); *\/ */
/*   vorticity_minus = Solution_Block.W[i][j][k-1].vorticity(Solution_Block.dWdx[i][j][k-1], */
/* 							  Solution_Block.dWdy[i][j][k-1], */
/* 							  Solution_Block.dWdz[i][j][k-1]); */
/*   DZ = Solution_Block.Grid.Cell[i][j][k].Xc.z - Solution_Block.Grid.Cell[i][j][k-1].Xc.z; */
/*   dvorx_dz_minus = (vorticity_plus.x - vorticity_minus.x)/DZ; */
/*   dvory_dz_minus = (vorticity_plus.y - vorticity_minus.y)/DZ; */
/*   dvorz_dz_minus = (vorticity_plus.z - vorticity_minus.z)/DZ; */


/*   //-----  Second derivatives -------// */

/*   DX = Solution_Block.Grid.Cell[i+1][j][k].Xc.x - Solution_Block.Grid.Cell[i-1][j][k].Xc.x; */
/*   d2vorx_dxx = (dvorx_dx_plus - dvorx_dx_minus)/DX; */
/*   d2vory_dxx = (dvory_dx_plus - dvory_dx_minus)/DX; */
/*   d2vorz_dxx = (dvorz_dx_plus - dvorz_dx_minus)/DX; */
    
/*   DY = Solution_Block.Grid.Cell[i][j+1][k].Xc.y - Solution_Block.Grid.Cell[i][j-1][k].Xc.y; */
/*   d2vorx_dyy = (dvorx_dy_plus - dvorx_dy_minus)/DY; */
/*   d2vory_dyy = (dvory_dy_plus - dvory_dy_minus)/DY; */
/*   d2vorz_dyy = (dvorz_dy_plus - dvorz_dy_minus)/DY; */

/*   DZ = Solution_Block.Grid.Cell[i][j][k+1].Xc.z - Solution_Block.Grid.Cell[i][j][k-1].Xc.z; */
/*   d2vorx_dzz = (dvorx_dz_plus - dvorx_dz_minus)/DZ; */
/*   d2vory_dzz = (dvory_dz_plus - dvory_dz_minus)/DZ; */
/*   d2vorz_dzz = (dvorz_dz_plus - dvorz_dz_minus)/DZ; */



  ///////////////////////////////////////////////////////////////////////////////////////
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


#endif /* _LES3D_HEXA_INCLUDED  */

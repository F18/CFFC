
#ifndef _TURBULENCE_AVERAGING_INCLUDED 
#define _TURBULENCE_AVERAGING_INCLUDED


/* Include required CFFC header files. */

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED




//-------------------------------------------------//
//                                                 //
//              Volume averaging                   //
//                                                 //
//-------------------------------------------------//


template<typename HEXA_BLOCK>
void Time_Averaging_of_Velocity_Field(HEXA_BLOCK *Solution_Block,
				      AdaptiveBlock3D_List &LocalSolnBlockList,
				      double &u_average,
				      double &v_average,
				      double &w_average) {

  double local_vol, Volume(ZERO);
  Vector3D vel;  vel.zero();
  
  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       vel += Solution_Block[p].W[i][j][k].v * local_vol;
	       Volume += local_vol; 
	   } /* endfor */
	} /* endfor */
      } /* endfor*/
    } /* endif */
  } /* endfor */

  Volume = CFFC_Summation_MPI(Volume);
  vel.x = CFFC_Summation_MPI(vel.x);
  vel.y = CFFC_Summation_MPI(vel.y);
  vel.z = CFFC_Summation_MPI(vel.z);
  u_average = vel.x/Volume;
  v_average = vel.y/Volume;
  w_average = vel.z/Volume;
}


template<typename HEXA_BLOCK>
void Time_Averaging_of_Solution(HEXA_BLOCK *Solution_Block,
				AdaptiveBlock3D_List &LocalSolnBlockList,
				const double &u_average, 
				const double &v_average,
				const double &w_average,
				double &sqr_u) {

  double vis, u_ave, v_ave, w_ave, local_vol, total_vol(ZERO), vis_ave(ZERO);
  double u_p(ZERO), v_p(ZERO), w_p(ZERO), ens(ZERO), eps_w(ZERO), eps_ss(ZERO);
  
  u_ave = u_average;
  v_ave = v_average;  
  w_ave = w_average;  
    
  for (int p = 0; p < LocalSolnBlockList.Nblk; p++) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl; i <= Solution_Block[p].ICu; ++i) {
        for (int j  = Solution_Block[p].JCl; j <= Solution_Block[p].JCu; ++j) {
           for (int k  = Solution_Block[p].KCl; k <= Solution_Block[p].KCu; ++k) {
                local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
                total_vol += local_vol;
                u_p += sqr(Solution_Block[p].W[i][j][k].v.x - u_ave) * local_vol;
                v_p += sqr(Solution_Block[p].W[i][j][k].v.y - v_ave) * local_vol;
                w_p += sqr(Solution_Block[p].W[i][j][k].v.z - w_ave) * local_vol;
                vis = ( Solution_Block[p].W[i][j][k].mu() + 
			Solution_Block[p].W[i][j][k].mu_t(Solution_Block[p].dWdx[i][j][k],
							  Solution_Block[p].dWdy[i][j][k],
							  Solution_Block[p].dWdz[i][j][k],
							  Solution_Block[p].Flow_Type, 
							  local_vol)
			) / Solution_Block[p].W[i][j][k].rho;
                vis_ave += vis*local_vol;
                ens += Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                              Solution_Block[p].dWdy[i][j][k],
                                                              Solution_Block[p].dWdz[i][j][k]) * local_vol;
                eps_w += 2.0*vis* Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                                         Solution_Block[p].dWdy[i][j][k],
                                                                         Solution_Block[p].dWdz[i][j][k])*local_vol;
                eps_ss += 2.0*vis*(sqr(Solution_Block[p].W[i][j][k].abs_strain_rate(Solution_Block[p].dWdx[i][j][k],
                                                                                    Solution_Block[p].dWdy[i][j][k],
                                                                                    Solution_Block[p].dWdz[i][j][k]))/ 
                          2.0)*local_vol;	      
	   } /* endfor */
	} /* endfor */
      } /* endfor*/
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  vis_ave = CFFC_Summation_MPI(vis_ave);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  w_p = CFFC_Summation_MPI(w_p);
  ens = CFFC_Summation_MPI(ens);
  eps_w = CFFC_Summation_MPI(eps_w);
  eps_ss = CFFC_Summation_MPI(eps_ss);

  sqr_u = u_p/total_vol;
  vis_ave = vis_ave/total_vol;
  u_p = u_p/total_vol;
  v_p = v_p/total_vol;
  w_p = w_p/total_vol;
  ens = ens/total_vol;
  eps_w = eps_w/total_vol;
  eps_ss = eps_ss/total_vol;
  
  double u_rms = sqrt((u_p + v_p + w_p)/3.0);
  double Taylor_scale, Kolmogorov_scale, Re_Taylor, L11; 
  double l_1, l_2;

  Kolmogorov_scale = pow(pow(vis_ave, THREE)/eps_ss, 0.25);

  if (ens == ZERO) {
    Taylor_scale = ZERO;
  } else {
    // Homogeneous isotropic
    Taylor_scale = sqrt(15.0*vis_ave*u_rms*u_rms/eps_ss);
    //cout << "\n =======> Taylor scale based on TKE dissipation: " << Taylor_scale;
    //cout << "\n =======> Taylor scale based on enstrophy: " << Taylor_scale;
  }

  Re_Taylor = u_rms*Taylor_scale/vis_ave;
  L11= 0.09*pow(0.5*u_rms*u_rms, 1.5)/eps_ss;

  if (eps_w > 0.0) {
    l_1 = 0.42*pow(u_rms, 3.0)/eps_w;
  } else {
    l_1 = 0.0;
  }

  if (eps_ss > 0.0) {
    l_2 = 0.42*pow(u_rms, 3.0)/eps_ss;
  } else {
    l_2 = 0.0;
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n ==========================================================================\n"; 
    cout << " Turbulent Statistics of Resolved Velocity Field (in Physical Space):\n";
    cout << "\n <u> = " << u_ave <<"  "<< "\t<v> = " << v_ave <<"  "<< "\t<w> = " << w_ave;
    cout << "\n eps_w = "<< eps_w <<"  "<< "\teps_ss = "<< eps_ss<<"  "<< "\tenstrophy = "<< ens;
    cout << "\n l_1 = "<< l_1 <<"  "<< "\tl_2 = " << l_2;
    cout << "\n nu = "<< vis_ave <<"  "<< "\tRe_Taylor = " << Re_Taylor << endl;
    cout << "\n ===> u_rms            = " << u_rms;
    cout << "\n ===> Taylor scale     = " << Taylor_scale;
    cout << "\n ===> Kolmogorov scale = " << Kolmogorov_scale;
    cout << "\n ===> L11              = " << L11 <<  endl;
    cout << " ==========================================================================" << endl;
  } /* endif */

}


// Total turbulence kinetic energy
template<typename HEXA_BLOCK>
double Total_TKE(HEXA_BLOCK *Solution_Block,
                 AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), u_p(ZERO), v_p(ZERO), w_p(ZERO);
  double u_ave, v_ave, w_ave, u_rms;

  Time_Averaging_of_Velocity_Field(Solution_Block, LocalSolnBlockList, u_ave, v_ave, w_ave);

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       u_p += sqr(Solution_Block[p].W[i][j][k].v.x - u_ave) * local_vol;
	       v_p += sqr(Solution_Block[p].W[i][j][k].v.y - v_ave) * local_vol;
	       w_p += sqr(Solution_Block[p].W[i][j][k].v.z - w_ave) * local_vol;
	  } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */
  
  total_vol = CFFC_Summation_MPI(total_vol);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  w_p = CFFC_Summation_MPI(w_p);
      
  u_p = u_p/total_vol;
  v_p = v_p/total_vol;
  w_p = w_p/total_vol;

  u_rms = sqrt((u_p + v_p + w_p)/3.0);
  
  // In 3D: k =  3*sqr(u_rms)/2
  return (3.0*u_rms*u_rms/2.0); 
}


// Total enstrophy
template <typename HEXA_BLOCK>
double Total_Enstrophy(HEXA_BLOCK *Solution_Block,
                       AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), ens(ZERO);

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       ens += Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                          Solution_Block[p].dWdy[i][j][k],
                                                          Solution_Block[p].dWdz[i][j][k]) * local_vol;
	   } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  ens = CFFC_Summation_MPI(ens);
 
  ens = ens/total_vol;
  
  return ens;
}


// Root mean square velocity (turbulence intensity)
template <typename HEXA_BLOCK>
double u_rms(HEXA_BLOCK *Solution_Block,
             AdaptiveBlock3D_List &LocalSolnBlockList) {

  double TKE;
  TKE = Total_TKE(Solution_Block, LocalSolnBlockList);

  return sqrt(2.0*TKE/3.0);
}


// Taylor scale of turbulence
template <typename HEXA_BLOCK>
double Taylor_Scale(HEXA_BLOCK *Solution_Block,
                    AdaptiveBlock3D_List &LocalSolnBlockList) {

  double enstrophy, taylor_scale, u_prime;
  enstrophy = Total_Enstrophy(Solution_Block, LocalSolnBlockList);

  u_prime = u_rms(Solution_Block, LocalSolnBlockList);
 
  if (enstrophy == ZERO) {
    taylor_scale = ZERO;
  } else {
    taylor_scale = sqrt(TWO*u_prime*u_prime/(TWO*enstrophy));
  }

  return taylor_scale;
}


// Volume-averaged kinematic viscosity
template <typename HEXA_BLOCK>
double Average_viscosity(HEXA_BLOCK *Solution_Block,
                         AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), vis(ZERO);

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       vis += ( Solution_Block[p].W[i][j][k].mu() + 
			Solution_Block[p].W[i][j][k].mu_t(Solution_Block[p].dWdx[i][j][k],
							  Solution_Block[p].dWdy[i][j][k],
							  Solution_Block[p].dWdz[i][j][k],
							  Solution_Block[p].Flow_Type, 
							  local_vol)
		      ) / Solution_Block[p].W[i][j][k].rho;
	   } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  vis = CFFC_Summation_MPI(vis);
    
  vis = vis/total_vol;

  return vis;
}




//------------------------------------------------------//
//                                                      //
//   Conditional volume averaging for premixed flames   //
//                                                      //
//------------------------------------------------------//

template<typename HEXA_BLOCK>
void Conditional_Averaging_of_Velocity_Field(HEXA_BLOCK *Solution_Block,
					     AdaptiveBlock3D_List &LocalSolnBlockList,
					     double &u_average,
					     double &v_average,
					     double &w_average) {

  double local_vol, Volume(ZERO), Yfuel_conditional(ZERO);
  Vector3D vel;  vel.zero();

  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;//Fresh_Fuel_Mass_Fraction;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       vel += Solution_Block[p].W[i][j][k].v * local_vol;
	       Volume += local_vol; 
	     } /* endif */
	   } /* endfor */
	} /* endfor */
      } /* endfor*/
    } /* endif */
  } /* endfor */

  Volume = CFFC_Summation_MPI(Volume);
  vel.x = CFFC_Summation_MPI(vel.x);
  vel.y = CFFC_Summation_MPI(vel.y);
  vel.z = CFFC_Summation_MPI(vel.z);
  u_average = vel.x/Volume;
  v_average = vel.y/Volume;
  w_average = vel.z/Volume;
}


template<typename HEXA_BLOCK>
void Conditional_Averaging_of_Solution(HEXA_BLOCK *Solution_Block,
				       AdaptiveBlock3D_List &LocalSolnBlockList,
				       const double &u_average, 
				       const double &v_average,
				       const double &w_average,
				       double &sqr_u) {

  double vis, u_ave, v_ave, w_ave, local_vol, total_vol(ZERO), vis_ave(ZERO);
  double u_p(ZERO), v_p(ZERO), w_p(ZERO), ens(ZERO), eps_w(ZERO), eps_ss(ZERO);
  double Yfuel_conditional(ZERO);
  
  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;//Fresh_Fuel_Mass_Fraction;
  u_ave = u_average;
  v_ave = v_average;  
  w_ave = w_average;  
    
  for (int p = 0; p < LocalSolnBlockList.Nblk; p++) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl; i <= Solution_Block[p].ICu; ++i) {
        for (int j  = Solution_Block[p].JCl; j <= Solution_Block[p].JCu; ++j) {
           for (int k  = Solution_Block[p].KCl; k <= Solution_Block[p].KCu; ++k) {
	      if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
                local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
                total_vol += local_vol;
                u_p += sqr(Solution_Block[p].W[i][j][k].v.x - u_ave) * local_vol;
                v_p += sqr(Solution_Block[p].W[i][j][k].v.y - v_ave) * local_vol;
                w_p += sqr(Solution_Block[p].W[i][j][k].v.z - w_ave) * local_vol;
		vis = ( Solution_Block[p].W[i][j][k].mu() + 
			Solution_Block[p].W[i][j][k].mu_t(Solution_Block[p].dWdx[i][j][k],
							  Solution_Block[p].dWdy[i][j][k],
							  Solution_Block[p].dWdz[i][j][k],
							  Solution_Block[p].Flow_Type, 
							  local_vol)
			) / Solution_Block[p].W[i][j][k].rho;
                vis_ave += vis*local_vol;
                ens += Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                              Solution_Block[p].dWdy[i][j][k],
                                                              Solution_Block[p].dWdz[i][j][k]) * local_vol;
                eps_w += 2.0*vis* Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                                         Solution_Block[p].dWdy[i][j][k],
                                                                         Solution_Block[p].dWdz[i][j][k])*local_vol;
                eps_ss += 2.0*vis*(sqr(Solution_Block[p].W[i][j][k].abs_strain_rate(Solution_Block[p].dWdx[i][j][k],
                                                                                    Solution_Block[p].dWdy[i][j][k],
                                                                                    Solution_Block[p].dWdz[i][j][k]))/ 
                          2.0)*local_vol;
	      } /* endif */
	   } /* endfor */
	} /* endfor */
      } /* endfor*/
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  vis_ave = CFFC_Summation_MPI(vis_ave);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  w_p = CFFC_Summation_MPI(w_p);
  ens = CFFC_Summation_MPI(ens);
  eps_w = CFFC_Summation_MPI(eps_w);
  eps_ss = CFFC_Summation_MPI(eps_ss);

  sqr_u = u_p/total_vol;
  vis_ave = vis_ave/total_vol;
  u_p = u_p/total_vol;
  v_p = v_p/total_vol;
  w_p = w_p/total_vol;
  ens = ens/total_vol;
  eps_w = eps_w/total_vol;
  eps_ss = eps_ss/total_vol;
  
  double u_rms = sqrt((u_p + v_p + w_p)/3.0);
  double Taylor_scale, Kolmogorov_scale, Re_Taylor, L11; 
  double l_1, l_2;

  Kolmogorov_scale = pow(pow(vis_ave, THREE)/eps_ss, 0.25);

  if (ens == ZERO) {
    Taylor_scale = ZERO;
  } else {
    // Homogeneous isotropic
    Taylor_scale = sqrt(15.0*vis_ave*u_rms*u_rms/eps_ss);
    //cout << "\n =======> Taylor scale based on TKE dissipation: " << Taylor_scale;
    //cout << "\n =======> Taylor scale based on enstrophy: " << Taylor_scale;
  }

  Re_Taylor = u_rms*Taylor_scale/vis_ave;
  L11= 0.09*pow(0.5*u_rms*u_rms, 1.5)/eps_ss;

  if (eps_w > 0.0) {
    l_1 = 0.42*pow(u_rms, 3.0)/eps_w;
  } else {
    l_1 = 0.0;
  }

  if (eps_ss > 0.0) {
    l_2 = 0.42*pow(u_rms, 3.0)/eps_ss;
  } else {
    l_2 = 0.0;
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n ==========================================================================\n"; 
    cout << " Turbulent Statistics of Resolved Velocity Field (in Physical Space):\n";
    cout << "\n <u> = " << u_ave <<"  "<< "\t<v> = " << v_ave <<"  "<< "\t<w> = " << w_ave;
    cout << "\n eps_w = "<< eps_w <<"  "<< "\teps_ss = "<< eps_ss<<"  "<< "\tenstrophy = "<< ens;
    cout << "\n l_1 = "<< l_1 <<"  "<< "\tl_2 = " << l_2;
    cout << "\n nu = "<< vis_ave <<"  "<< "\tRe_Taylor = " << Re_Taylor << endl;
    cout << "\n ===> u_rms            = " << u_rms;
    cout << "\n ===> Taylor scale     = " << Taylor_scale;
    cout << "\n ===> Kolmogorov scale = " << Kolmogorov_scale;
    cout << "\n ===> L11              = " << L11 <<  endl;
    cout << " ==========================================================================" << endl;
  } /* endif */

}


// Conditional total turbulence kinetic energy
template<typename HEXA_BLOCK>
double Conditional_Total_TKE(HEXA_BLOCK *Solution_Block,
			     AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), u_p(ZERO), v_p(ZERO), w_p(ZERO);
  double u_ave, v_ave, w_ave, u_rms;
  double Yfuel_conditional(ZERO); 

  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;

  Conditional_Averaging_of_Velocity_Field(Solution_Block, LocalSolnBlockList, u_ave, v_ave, w_ave);

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       u_p += sqr(Solution_Block[p].W[i][j][k].v.x - u_ave) * local_vol;
	       v_p += sqr(Solution_Block[p].W[i][j][k].v.y - v_ave) * local_vol;
	       w_p += sqr(Solution_Block[p].W[i][j][k].v.z - w_ave) * local_vol;
	   } /* endif */
	  } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */
  
  total_vol = CFFC_Summation_MPI(total_vol);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  w_p = CFFC_Summation_MPI(w_p);
      
  u_p = u_p/total_vol;
  v_p = v_p/total_vol;
  w_p = w_p/total_vol;

  u_rms = sqrt((u_p + v_p + w_p)/3.0);
  
  // In 3D: k =  3*sqr(u_rms)/2
  return (3.0*u_rms*u_rms/2.0); 
}


// Conditional total enstrophy
template <typename HEXA_BLOCK>
double Conditional_Total_Enstrophy(HEXA_BLOCK *Solution_Block,
				   AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), ens(ZERO);
  double Yfuel_conditional(ZERO);
   
  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       ens += Solution_Block[p].W[i][j][k].Enstrophy(Solution_Block[p].dWdx[i][j][k],
                                                          Solution_Block[p].dWdy[i][j][k],
                                                          Solution_Block[p].dWdz[i][j][k]) * local_vol;
	   } /* endif */
	  } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  ens = CFFC_Summation_MPI(ens);
 
  ens = ens/total_vol;
  
  return ens;
}


// Conditional root mean square velocity (turbulence intensity)
template <typename HEXA_BLOCK>
double Conditional_u_rms(HEXA_BLOCK *Solution_Block,
			 AdaptiveBlock3D_List &LocalSolnBlockList) {

  double TKE;
  TKE = Conditional_Total_TKE(Solution_Block, LocalSolnBlockList);

  return sqrt(2.0*TKE/3.0);
}


// Conditional Taylor scale of turbulence
template <typename HEXA_BLOCK>
double Conditional_Taylor_Scale(HEXA_BLOCK *Solution_Block,
				AdaptiveBlock3D_List &LocalSolnBlockList) {

  double enstrophy, taylor_scale, u_prime;
  enstrophy = Conditional_Total_Enstrophy(Solution_Block, LocalSolnBlockList);

  u_prime = Conditional_u_rms(Solution_Block, LocalSolnBlockList);
 
  if (enstrophy == ZERO) {
    taylor_scale = ZERO;
  } else {
    taylor_scale = sqrt(TWO*u_prime*u_prime/(TWO*enstrophy));
  }

  return taylor_scale;
}


// Conditional volume-averaged kinematic viscosity
template <typename HEXA_BLOCK>
double Conditional_Average_viscosity(HEXA_BLOCK *Solution_Block,
				     AdaptiveBlock3D_List &LocalSolnBlockList) {

  double local_vol, total_vol(ZERO), vis(ZERO);
  double Yfuel_conditional(ZERO);
  
  //Conditional average on fresh gas
  Yfuel_conditional = 0.95*0.05518;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
           for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	     if (Solution_Block[p].W[i][j][k].spec[0].c >= Yfuel_conditional) {
	       local_vol = Solution_Block[p].Grid.Cell[i][j][k].V;
	       total_vol += local_vol;
	       vis += ( Solution_Block[p].W[i][j][k].mu() + 
			Solution_Block[p].W[i][j][k].mu_t(Solution_Block[p].dWdx[i][j][k],
							  Solution_Block[p].dWdy[i][j][k],
							  Solution_Block[p].dWdz[i][j][k],
							  Solution_Block[p].Flow_Type, 
							  local_vol)
		      ) / Solution_Block[p].W[i][j][k].rho;
	     } /* endif */
	   } /* endfor */
	} /* endfor */
      } /* endfor */
    } /* endif */
  } /* endfor */

  total_vol = CFFC_Summation_MPI(total_vol);
  vis = CFFC_Summation_MPI(vis);
    
  vis = vis/total_vol;

  return vis;
}

// Maximum and minimum cell volumes in the domain
template<typename HEXA_BLOCK>
void Max_and_Min_Cell_Volumes(HEXA_BLOCK *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList) {

  double min_volume(1E20), max_volume(0.0), cell_volume;

  for (int p = 0 ; p <= LocalSolnBlockList.Nblk-1 ; p++ ) {
    if (LocalSolnBlockList.Block[p].used == ADAPTIVEBLOCK3D_USED) {
      for (int i = Solution_Block[p].ICl ; i <= Solution_Block[p].ICu ; i++) {
        for (int j = Solution_Block[p].JCl ; j <= Solution_Block[p].JCu ; j++) {
	  for (int k = Solution_Block[p].KCl ; k <= Solution_Block[p].KCu ; k++) {
	    cell_volume = Solution_Block[p].Grid.Cell[i][j][k].V;
	    if ( cell_volume > max_volume ) max_volume = cell_volume;
	    if ( cell_volume < min_volume ) min_volume = cell_volume;
	  }
	}
      }
    }
  }

  max_volume = CFFC_Maximum_MPI(max_volume);
  min_volume = CFFC_Minimum_MPI(min_volume);

  if ( CFFC_Primary_MPI_Processor() ) {
    cout << "\n -------------------------------------------"; 
    cout << "\n ===> Maximum cell volume  = " << max_volume;
    cout << "\n ===> Minimum cell volume  = " << min_volume;
    cout << "\n ===> Maximum cell size    = " << pow(max_volume, 1.0/3.0);
    cout << "\n ===> Minimum cell size    = " << pow(min_volume, 1.0/3.0);
    cout << "\n -------------------------------------------\n";
  }

}

/********************************************************
 *       Burning rate                                   *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double Turbulent_Burning_Rate(Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> *Solution_Block,
			      AdaptiveBlock3D_List &LocalSolnBlockList,
			      Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {

  cout << "\n Explicit Specialization of "
       << "Turbulent_Burning_Rate() "
       << "required.\n";

  return ZERO;
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
 *  Return the magnitude of vorticity at specified node *
 ********************************************************/
template<typename HEXA_BLOCK>
double vorticity_n(HEXA_BLOCK &Soln_Blk,
		   const int &i, 
		   const int &j, 
		   const int &k){

  return Trilinear_Interpolation(
    Soln_Blk.Grid.Cell[i-1][j][k].Xc, Soln_Blk.W[i-1][j][k].vorticity(Soln_Blk.dWdx[i-1][j][k],
								      Soln_Blk.dWdy[i-1][j][k],
								      Soln_Blk.dWdz[i-1][j][k]).abs(),
    Soln_Blk.Grid.Cell[i][j][k].Xc, Soln_Blk.W[i][j][k].vorticity(Soln_Blk.dWdx[i][j][k],
								  Soln_Blk.dWdy[i][j][k],
								  Soln_Blk.dWdz[i][j][k]).abs(),
    Soln_Blk.Grid.Cell[i][j-1][k].Xc, Soln_Blk.W[i][j-1][k].vorticity(Soln_Blk.dWdx[i][j-1][k],
								      Soln_Blk.dWdy[i][j-1][k],
								      Soln_Blk.dWdz[i][j-1][k]).abs(),
    Soln_Blk.Grid.Cell[i-1][j-1][k].Xc, Soln_Blk.W[i-1][j-1][k].vorticity(Soln_Blk.dWdx[i-1][j-1][k],
									  Soln_Blk.dWdy[i-1][j-1][k],
									  Soln_Blk.dWdz[i-1][j-1][k]).abs(),
    Soln_Blk.Grid.Cell[i-1][j][k-1].Xc, Soln_Blk.W[i-1][j][k-1].vorticity(Soln_Blk.dWdx[i-1][j][k-1],
									  Soln_Blk.dWdy[i-1][j][k-1],
									  Soln_Blk.dWdz[i-1][j][k-1]).abs(),
    Soln_Blk.Grid.Cell[i][j][k-1].Xc, Soln_Blk.W[i][j][k-1].vorticity(Soln_Blk.dWdx[i][j][k-1],
								      Soln_Blk.dWdy[i][j][k-1],
								      Soln_Blk.dWdz[i][j][k-1]).abs(),
    Soln_Blk.Grid.Cell[i][j-1][k-1].Xc, Soln_Blk.W[i][j-1][k-1].vorticity(Soln_Blk.dWdx[i][j-1][k-1],
									  Soln_Blk.dWdy[i][j-1][k-1],
									  Soln_Blk.dWdz[i][j-1][k-1]).abs(),
    Soln_Blk.Grid.Cell[i-1][j-1][k-1].Xc, Soln_Blk.W[i-1][j-1][k-1].vorticity(Soln_Blk.dWdx[i-1][j-1][k-1],
									      Soln_Blk.dWdy[i-1][j-1][k-1],
									      Soln_Blk.dWdz[i-1][j-1][k-1]).abs(),
    Soln_Blk.Grid.Node[i][j][k].X);

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

#endif // _TURBULENCE_AVERAGING_INCLUDED

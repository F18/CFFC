#ifndef _LESPREMIXED2D_TURB_INCLUDED
#include "LESPremixed2DTurbInit.h"
#endif

#ifndef _LESPREMIXED2D_QUAD_INCLUDED
#include "LESPremixed2DQuad.h"
#endif



/*======================================*\
        Single block functions
\*======================================*/  


// Calculate the area of the block
double Total_Block_Area(LESPremixed2D_Quad_Block &SolnBlk) {
  double total_A = ZERO;

  for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
    for (int j  = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      total_A += SolnBlk.Grid.Cell[i][j].A;      
    }
  }

  return total_A;
}




/*======================================*\
      Multiple block functions          
\*======================================*/


// Total area
double Total_Area(LESPremixed2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List) {

  double Area = ZERO;

  for (int i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      Area += Total_Block_Area(Soln_ptr[i]);
    }
  }

  Area = CFFC_Summation_MPI(Area);

  return Area;
} 
 
              
// Maximum wave number
double Max_k(LESPremixed2D_Input_Parameters &Input_Parameters) {
  double k_max, abs_k, k1, k2;
  int Nx = Input_Parameters.Number_of_Cells_Idir;
  int Ny = Input_Parameters.Number_of_Cells_Jdir;

  k_max = ZERO;

  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      if (abs_k > k_max) {
        k_max = abs_k;
      }
    }
  }
  return k_max;
} 


// Obtain the area-averaged velocity
void Average_V(LESPremixed2D_Quad_Block *Soln_ptr,
               AdaptiveBlock2D_List &Soln_Block_List,
               LESPremixed2D_Input_Parameters &Input_Parameters,
               double &u_average, double &v_average) {

  double local_A, Area = ZERO;
  double Yfuel_conditional = ZERO;
  Vector2D vel;  vel.zero();
  
  //Conditional average on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }

  for (int p = 0 ; p <= Soln_Block_List.Nblk-1 ; p++ ) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; i++) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; j++) {
          if (Soln_ptr[p].W[i][j].spec[0].c >= Yfuel_conditional) {
	    local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	    vel += Soln_ptr[p].W[i][j].v * local_A;
	  }
        }
      }
      Area += Total_Block_Area(Soln_ptr[p]);  
    }
  }

  Area = CFFC_Summation_MPI(Area);
  vel.x = CFFC_Summation_MPI(vel.x);
  vel.y = CFFC_Summation_MPI(vel.y);
    
  u_average = vel.x/Area;
  v_average = vel.y/Area;
  //cout << "<u> = "<<u_average<<"  <v> = " <<v_average<<flush<<endl;
  
}



// Obtain the area-averaged kinematic viscosity
double Average_viscosity(LESPremixed2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_A, Area = ZERO, vis = ZERO;
  double Yfuel_conditional = ZERO;
  
   //Conditional average on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }

  for (int p = 0 ; p <= Soln_Block_List.Nblk-1 ; p++ ) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; i++) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; j++) {
          if (Soln_ptr[p].W[i][j].spec[0].c >= Yfuel_conditional) {
	    local_A = Soln_ptr[p].Grid.Cell[i][j].A;
            // Rescale viscosity dividing by the thickening factor and the wrinkling
            // factor,i.e. calculate the real viscosity.
#ifdef THICKENED_FLAME_ON
	    vis += Soln_ptr[p].W[i][j].mu()*local_A/(Soln_ptr[p].W[i][j].rho*Soln_ptr[p].W[i][j].flame.TF*Soln_ptr[p].W[i][j].flame.WF);
#else
	    vis += Soln_ptr[p].W[i][j].mu()*local_A/Soln_ptr[p].W[i][j].rho;
#endif
	  }
        }
      }
      Area += Total_Block_Area(Soln_ptr[p]);  
    }
  }

  Area = CFFC_Summation_MPI(Area);
  vis = CFFC_Summation_MPI(vis);
    
  vis = vis/Area;

  return vis;
}



// Total turbulence kinetic energy
double Total_TKE(LESPremixed2D_Quad_Block *Soln_ptr,
                 AdaptiveBlock2D_List &Soln_Block_List,
                 LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_A, total_A = ZERO, u_p = ZERO, v_p = ZERO;
  double u_ave, v_ave, u_rms;
  double Yfuel_conditional = ZERO; 

   //Conditional average on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }

  Average_V(Soln_ptr, Soln_Block_List, Input_Parameters, u_ave, v_ave);

  for (int p = 0; p <= Soln_Block_List.Nblk-1; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          if (Soln_ptr[p].W[i][j].spec[0].c >= Yfuel_conditional) {
	    local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	    total_A += local_A;
	    u_p += sqr(Soln_ptr[p].W[i][j].v.x - u_ave) * local_A;
	    v_p += sqr(Soln_ptr[p].W[i][j].v.y - v_ave) * local_A;
	  }
        }
      }
    }
  }
  
  total_A = CFFC_Summation_MPI(total_A);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
      
  u_p = u_p/total_A;
  v_p = v_p/total_A;
  u_rms = sqrt(0.5*(u_p + v_p));
  
  // In 2D: k = 2 sqr(u_rms)/2
  return (u_rms*u_rms); 
}



// Total enstrophy
double Total_Enstrophy(LESPremixed2D_Quad_Block *Soln_ptr,
                       AdaptiveBlock2D_List &Soln_Block_List,
                       LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_A, total_A = ZERO, ens = ZERO;
  double Yfuel_conditional = ZERO;
   
  //Conditional average on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }

  for (int p = 0; p <= Soln_Block_List.Nblk-1; p++) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {  
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          if (Soln_ptr[p].W[i][j].spec[0].c >= Yfuel_conditional) {
	    local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	    total_A += local_A;
	    ens += Soln_ptr[p].enstrophy(i, j) * local_A;      
	  }
        }
      }
    }
  }

  total_A = CFFC_Summation_MPI(total_A);
  ens = CFFC_Summation_MPI(ens);
 
  //CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  ens = ens/total_A;
  
  return ens;
}



// Root mean square velocity (turbulence intensity)
double u_rms(LESPremixed2D_Quad_Block *Soln_ptr,
	     AdaptiveBlock2D_List &Soln_Block_List,
	     LESPremixed2D_Input_Parameters &Input_Parameters) {

  double TKE;
  TKE = Total_TKE(Soln_ptr, Soln_Block_List, Input_Parameters);

  return sqrt(TKE); // 2D
}


// Taylor scale of turbulence
double Taylor_Scale(LESPremixed2D_Quad_Block *Soln_ptr,
		    AdaptiveBlock2D_List &Soln_Block_List,
		    LESPremixed2D_Input_Parameters &Input_Parameters) {

  double enstrophy, taylor_scale, u_prime;
  enstrophy = Total_Enstrophy(Soln_ptr, Soln_Block_List, Input_Parameters);

  u_prime = u_rms(Soln_ptr, Soln_Block_List, Input_Parameters);
 
  if (enstrophy == ZERO) {
    taylor_scale = ZERO;
  } else {
    taylor_scale = sqrt(TWO*u_prime*u_prime/(TWO*enstrophy));
  }

  return taylor_scale;
}



/*----------------------------------------------------------------*\
   Turbulent burning rate based on fuel and simplified for a 
   rectangular box.
\*----------------------------------------------------------------*/
double Turbulent_Burning_Rate(LESPremixed2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_A, Yf_u, rho_u, Ly, burning_rate = ZERO;
  Yf_u = Input_Parameters.Fresh_Fuel_Mass_Fraction;
  rho_u = Input_Parameters.Fresh_Density;
  Ly = Input_Parameters.Box_Height; 

  for (int p = 0; p <= Soln_Block_List.Nblk-1; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {  
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	  if (Soln_ptr[p].Flow_Type == FLOWTYPE_LAMINAR_C_FSD ||
              Soln_ptr[p].Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD ||
              Soln_ptr[p].Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
              Soln_ptr[p].Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY ||
              Soln_ptr[p].Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE ) {
	    burning_rate +=  Soln_ptr[p].W[i][j].scalar[1]*local_A*Soln_ptr[p].W[i][j].rho;
	    //	    	    burning_rate += Soln_ptr[p].W[i][j].Reaction_Rate_Fsd( Soln_ptr[p].dWdx[i][j],Soln_ptr[p].dWdy[i][j])*local_A;
	  }else{
          // Rate of consumption of fuel
	  burning_rate += Soln_ptr[p].W[i][j].Sw(Soln_ptr[p].W[i][j].React.reactset_flag,
                                                 Soln_ptr[p].Flow_Type).rhospec[0].c*local_A;      
        }
      }
    }
  }
  }
  burning_rate = CFFC_Summation_MPI(burning_rate);
  if (Input_Parameters.FlowType == FLOWTYPE_LAMINAR_C_FSD ||
      Input_Parameters.FlowType == FLOWTYPE_LAMINAR_NGT_C_FSD ||
      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY ||
      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE ) {
    burning_rate = burning_rate*Input_Parameters.laminar_flame_speed/Ly;//(rho_u*Ly);  //(rho_u*Yf_u*Ly);
  }else{
  burning_rate = -burning_rate/(rho_u*Yf_u*Ly);
  }  
  return burning_rate;
}

double Turbulent_Burning_Rate_Progvar(LESPremixed2D_Quad_Block *Soln_ptr,
			              AdaptiveBlock2D_List &Soln_Block_List,
			              LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_A, Yf_u, rho_u, Ly, lam_speed_fsd, turbulent_burning_rate_prog = ZERO;
  Yf_u = Input_Parameters.Fresh_Fuel_Mass_Fraction;
  rho_u = Input_Parameters.Fresh_Density;
  lam_speed_fsd = Input_Parameters.laminar_flame_speed;
  Ly = Input_Parameters.Box_Height; 

  for (int p = 0; p <= Soln_Block_List.Nblk-1; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {  
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	    Tensor2D strain_rate;
	    strain_rate = Soln_ptr[p].W[i][j].Strain_Rate(Soln_ptr[p].dWdx[i][j], Soln_ptr[p].dWdy[i][j], 
					  	          Soln_ptr[p].Flow_Type, Soln_ptr[p].Axisymmetric, 
						          Soln_ptr[p].Grid.Cell[i][j].Xc);  
 	  if (Soln_ptr[p].Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC ||
              Soln_ptr[p].Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC ) {
          turbulent_burning_rate_prog += Soln_ptr[p].W[i][j].Reaction_Rate_Algebraic( Soln_ptr[p].dWdx[i][j],Soln_ptr[p].dWdy[i][j],Soln_ptr[p].Flow_Type,strain_rate)*local_A;
 	  }else{
	  turbulent_burning_rate_prog += Soln_ptr[p].W[i][j].Reaction_Rate_Progvar( Soln_ptr[p].dWdx[i][j],Soln_ptr[p].dWdy[i][j])*local_A;//Soln_ptr[p].W[i][j].Reaction_Rate_Fsd_Algebraic(Soln_ptr[p].dWdx[i][j],Soln_ptr[p].dWdy[i][j],local_A,Soln_ptr[p].Flow_Type)*local_A;//Soln_ptr[p].W[i][j].Reaction_Rate_Progvar( Soln_ptr[p].dWdx[i][j],Soln_ptr[p].dWdy[i][j])*local_A;
	  }
	 }
        }
      }
    }
  
  turbulent_burning_rate_prog = CFFC_Summation_MPI(turbulent_burning_rate_prog);
  turbulent_burning_rate_prog = turbulent_burning_rate_prog/(rho_u*Ly);//*lam_speed_fsd/Ly;//(rho_u*Ly);  //(rho_u*Yf_u*Ly);
  return turbulent_burning_rate_prog;
}

/*----------------------------------------------------------------*\
   Total mass of species k per unit transversal length 
\*----------------------------------------------------------------*/
double Total_Species_Mass(LESPremixed2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  const int &k) {

  double local_A, species_mass = ZERO;

  for (int p = 0; p <= Soln_Block_List.Nblk-1; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {  
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          local_A = Soln_ptr[p].Grid.Cell[i][j].A;
	  // rho * Y[k] * A
	  species_mass += Soln_ptr[p].W[i][j].rho*Soln_ptr[p].W[i][j].spec[k].c*local_A;   
        }
      }
    }
  }

  species_mass = CFFC_Summation_MPI(species_mass);
  if (CFFC_Primary_MPI_Processor()) {
    cout  << "\n Mass of species[" << k << "] = " << species_mass;
  }
  return species_mass;
}


// Average x position of the flame
double Average_Flame_Position(LESPremixed2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LESPremixed2D_Input_Parameters &Input_Parameters) {

  double local_progress_variable, Yf_u, Yf_b, Yf;
  double local_A, Area, mean_x;
  
  Yf_u = Input_Parameters.Fresh_Fuel_Mass_Fraction;
  Yf_b = Input_Parameters.Burnt_Fuel_Mass_Fraction;
  
  Area = ZERO;
  mean_x = ZERO;
  for (int p = 0; p <= Soln_Block_List.Nblk-1; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {  
      for (int i = Soln_ptr[p].ICl ; i <= Soln_ptr[p].ICu ; ++i) {
        for (int j = Soln_ptr[p].JCl ; j <= Soln_ptr[p].JCu ; ++j) {
          Yf = Soln_ptr[p].W[i][j].spec[0].c;
          local_progress_variable = (Yf - Yf_u)/(Yf_b - Yf_u);
          if (local_progress_variable > 0.03  &&  local_progress_variable < 0.97) {
            local_A = Soln_ptr[p].Grid.Cell[i][j].A;
            Area += local_A;
            mean_x += Soln_ptr[p].Grid.Cell[i][j].Xc.x*local_A;
	  }
        }
      }
    }
  }

  Area = CFFC_Summation_MPI(Area);
  mean_x = CFFC_Summation_MPI(mean_x);

  mean_x = mean_x/Area;

  return mean_x;
}



// Turbulence statistics in physical space
void Average(LESPremixed2D_Quad_Block *Soln_ptr,
             AdaptiveBlock2D_List &Soln_Block_List,
             LESPremixed2D_Input_Parameters &Input_Parameters,
             const double &u_average, const double &v_average,
             double &sqr_u) {

  double vis, u_ave, v_ave, local_A, total_A = ZERO, vis_ave = ZERO;
  double u_p=ZERO, v_p=ZERO, ens=ZERO, eps_w=ZERO;  //eps_ss=ZERO;
  double Yfuel_conditional = ZERO;
  LESPremixed2D_pState dWdx, dWdy;
  
  //Conditional average on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }
      
  u_ave = u_average;
  v_ave = v_average;  
    
  for (int p = 0; p < Soln_Block_List.Nblk; p++) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {
      int q = p; 
      for (int i = Soln_ptr[q].ICl; i <= Soln_ptr[q].ICu; ++i) {
        for (int j  = Soln_ptr[q].JCl; j <= Soln_ptr[q].JCu; ++j) {
          if (Soln_ptr[q].W[i][j].spec[0].c >= Yfuel_conditional) {
            vis = Soln_ptr[q].W[i][j].mu()/(Soln_ptr[q].W[i][j].rho);
	    // Rescale viscosity dividing by the thickening factor and the wrinkling
            // factor,i.e. calculate the real viscosity.
#ifdef THICKENED_FLAME_ON
	    vis = vis/(Soln_ptr[q].W[i][j].flame.TF*Soln_ptr[q].W[i][j].flame.WF);
#endif
            local_A = Soln_ptr[q].Grid.Cell[i][j].A;
            total_A += local_A;
            vis_ave += vis*local_A;
            u_p += sqr(Soln_ptr[q].W[i][j].v.x - u_ave) * local_A;
            v_p += sqr(Soln_ptr[q].W[i][j].v.y - v_ave) * local_A;
            ens += Soln_ptr[q].enstrophy(i, j) * local_A;
            eps_w += 2.0*vis* Soln_ptr[q].enstrophy(i, j)*local_A;
            dWdx = Soln_ptr[q].dWdx[i][j];
            dWdy = Soln_ptr[q].dWdy[i][j];
            //eps_ss += 2.0*vis* (sqr(Soln_ptr[q].W[i][j].abs_strain_rate(dWdx, dWdy))/2.0)
            //        * local_A;
	  }
	}
      }
    }
  }
  
  total_A = CFFC_Summation_MPI(total_A);
  vis_ave = CFFC_Summation_MPI(vis_ave);
  u_p = CFFC_Summation_MPI(u_p);
  v_p = CFFC_Summation_MPI(v_p);
  ens = CFFC_Summation_MPI(ens);
  eps_w = CFFC_Summation_MPI(eps_w);
  //eps_ss = CFFC_Summation_MPI(eps_ss);

  sqr_u = u_p/total_A;

   
  vis_ave = vis_ave/total_A;
  u_p = u_p/total_A;
  v_p = v_p/total_A;
  ens = ens/total_A;
  eps_w = eps_w/total_A;
  //eps_ss = eps_ss/total_A;
  
  double u_rms = sqrt(0.5*(u_p + v_p));
  double Taylor_scale, Kolmogorov_scale, Re_Taylor; 
  double l_1, l_2;

  if (ens == ZERO) {
    Taylor_scale = ZERO;
  } else {
    Taylor_scale = sqrt(TWO*u_rms*u_rms/(TWO*ens));
  }

  Re_Taylor =  u_rms*Taylor_scale/vis_ave;
  Kolmogorov_scale = pow(pow(vis_ave, THREE)/eps_w, 0.25);
double L11= 0.09*pow(u_rms*u_rms, 1.5)/eps_w;

  if (eps_w > 0.0) {
    l_1 = 0.42*pow(u_rms, 3.0)/eps_w;
  } else {
    l_1 = 0.0;
  }
  // if (eps_ss > 0.0) {
//     l_2 = 0.42*pow(u_rms, 3.0)/eps_ss;
//   } else {
//     l_2 = 0.0;
//   }

  if (CFFC_Primary_MPI_Processor()) {
    // Output
    cout << "\n\n ==========================================================================\n"; 
    cout << " In physical space:\n";
    cout << "\n <u^2> = "<< u_p <<"  "<< "<v^2> = "<< v_p <<"  "
	 << "u_rms  = " << u_rms <<"  "
	 << "\n <u> = " << u_ave <<"  "<< "<v> = " << v_ave <<"  "
	 << "ens = "<< ens <<"  " 
	 << "\n eps_w = "<< eps_w  //<<"  "<< "eps_ss = "<< eps_ss <<"  "
	 << "l_1 = "<< l_1 <<"  "<< "l_2 = " << l_2 <<"  "
	 << "\n vis = "<< vis_ave << "  Re_Taylor = " << Re_Taylor <<"  "
	 << "\n Taylor_scale = " << Taylor_scale <<"  "<<"L11=  "<<L11<<" " 
	 << "Kolmogorov_scale = " << Kolmogorov_scale << endl;
    cout << " ==========================================================================" << endl;
  } 
  
}



// Velocity fluctuations for the initialization of the turbulent field
void Velocity_Fluctuations(Grid2D_Quad_Block   **InitMeshBlk,
			   QuadTreeBlock_DataStructure &QuadTree,
			   LESPremixed2D_Input_Parameters &Input_Parameters) {

  int Nx, Ny, ny;
  double Lx, Ly; 
  Lx = Input_Parameters.Box_Width;
  Ly = Input_Parameters.Box_Height;
  Nx = Input_Parameters.Number_of_Cells_Idir;
  Ny = Input_Parameters.Number_of_Cells_Jdir;
  ny = Ny/2+1;
  
  //double        scale = 1.0/(Nx*Ny);  // Scaling factor for the real to complex transform
  double          *u, *v;               // Arrays to store the velocity fluctuations in physical space
  fftw_complex    *uu, *vv;             // Arrays to store the velocity fluctuations in Fourier space
  fftw_plan       physical;
  

  // Allocation of arrays used in the transforms
  u = (double *) malloc(Nx*Ny * sizeof(double));
  v = (double *) malloc(Nx*Ny * sizeof(double));
  uu = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
  vv = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
     
 
  int seed = 5;                    // fixed seed for the random number generator
  //int seed = time(NULL);         // assigns the current time to the seed
  srand48(seed);                   // changes the seed for drand48()
  int iconj;                       // position of the conjugate complex for the i index
  double k1, k2;                   // wave numbers
  double *EE_x, *EE_y, *Crf;       // spectra and rescaling coefficient  
  ofstream OutInitSpectrum;  


  // Allocation of componentes of the spectra
  double k_max = Max_k(Input_Parameters);
  cout <<"\n kmax = "<< k_max << endl;
  EE_x = new double[int(k_max)+1];
  EE_y = new double[int(k_max)+1];
  Crf  = new double[int(k_max)+1];

  for (int p=0; p<=int(k_max); p++) {
    EE_x[p] = 0.0;
    EE_y[p] = 0.0;
  }


  Complex aa;
  double abs_k, theta;
  for (int i=0; i<Nx; i++) {
    iconj = (i==0  ?  0 : Nx-i);
    for (int j=1; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      // Rogallo's function
      theta =  2.0*PI*random_double();   // Random number between 0 and 2*PI
      aa = alpha_Rogallo(abs_k, theta, Input_Parameters.i_Spectrum);

      // Assign spectral components of velocity
      uu[i*ny+j][0] = real(aa);
      uu[i*ny+j][1] = imag(aa);
      vv[i*ny+j][0] = uu[i*ny+j][0];
      vv[i*ny+j][1] = uu[i*ny+j][1];
                 
      if (i==0  &&  j==0) {
        uu[i*ny+j][0] = 0.0;
        uu[i*ny+j][1] = 0.0;
        vv[i*ny+j][0] = uu[i*ny+j][0];
	vv[i*ny+j][1] = uu[i*ny+j][1];
      } else {
        uu[i*ny+j][0] = k2*uu[i*ny+j][0]/abs_k;
        uu[i*ny+j][1] = k2*uu[i*ny+j][1]/abs_k;
        vv[i*ny+j][0] = -k1*vv[i*ny+j][0]/abs_k;
        vv[i*ny+j][1] = -k1*vv[i*ny+j][1]/abs_k;
      }
            
      // j=0;
      if (j==0 && i>Nx/2) {
        uu[i*ny][0] = uu[iconj*ny][0];
        uu[i*ny][1] = -uu[iconj*ny][1];
        vv[i*ny][0] = vv[iconj*ny][0];
        vv[i*ny][1] = -vv[iconj*ny][1];
      }

      // corners
      if (i==0   || i==Nx/2) {
        uu[i*ny][1] = 0.0;
        uu[i*ny+Ny/2][1] = 0.0;
        vv[i*ny][1] = 0.0;
        vv[i*ny+Ny/2][1] = 0.0;
      }
      

      // Compute the divergence of the fluctuations in Fourier space
      double div_re, div_im;
      div_re = k1*uu[i*ny+j][0] + k2*vv[i*ny+j][0];
      div_im = k1*uu[i*ny+j][1] + k2*vv[i*ny+j][1];
      if ( fabs(div_re) > 1.0E-6  ||  fabs(div_im) > 1.0E-6) {
        cout << "\n Warning: Divergence greater than tolerance ("
             << div_re << ", " << div_im << ")" << endl;
      }

     
      // Compute the initial spectra
      if (j == 0  ||  j == Ny/2) {
        EE_x[int(abs_k)] += 1.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 1.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      } else {
        EE_x[int(abs_k)] += 2.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 2.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      }
          
    }
    //cout << endl;
  }


  double int_x = 0.0, int_y = 0.0; // spectra checks (integration of the spectra)
  double ens = 0.0;                // enstrophy
  double l_11 = 0.0;               // L_11 integral length scale
  double EE, Eo, l;


  //------ Calculate the scaling factor of the velocity fluctuations  -------//
  for (int k=0; k<=int(k_max); k++) {
    EE = 0.5*(EE_x[k]+EE_y[k]);
    Eo = Energy_Spectrum(double(k), Input_Parameters.i_Spectrum);
    if (EE == ZERO) {
      Crf[k] = ZERO; 
    } else {
      Crf[k] = sqrt(Eo/EE);
    }
  }
  
  //-----------      Rescale fluctuations       -------------//
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      // Rescaling
      uu[i*ny+j][0] = Crf[int(abs_k)]*uu[i*ny+j][0];
      uu[i*ny+j][1] = Crf[int(abs_k)]*uu[i*ny+j][1]; 
      vv[i*ny+j][0] = Crf[int(abs_k)]*vv[i*ny+j][0];
      vv[i*ny+j][1] = Crf[int(abs_k)]*vv[i*ny+j][1];   
    }
  }

  //----------    Calculate the rescaled energy spectrum   ---------//
  // Reset EE_x and EE_y to ZERO
  for (int p=0; p<=int(k_max); p++) {
    EE_x[p] = 0.0;
    EE_y[p] = 0.0;
  }

  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      if (j == 0  ||  j == Ny/2) {
        EE_x[int(abs_k)] += 1.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 1.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      } else {
        EE_x[int(abs_k)] += 2.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 2.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      }
    }
  }

  
  // Output to Init_Spectrum file and checks
  OutInitSpectrum.open("Init_Spectrum.dat", ios::out);
  if(OutInitSpectrum.fail()){ 
    cerr<<"\nError opening file: Init_Spectrum.dat" << endl;
    exit(1); 
  }	  

  for (int k=0; k<=int(k_max); k++) {
    int_x += EE_x[k];
    int_y += EE_y[k];
    EE = 0.5*(EE_x[k]+EE_y[k]);
    if (k==0) {
      l_11 += 0.0; 
    } else {
      l_11 += EE/double(k);
    } 
    ens += EE*sqr(double(k));
    OutInitSpectrum << k <<"  "<< EE_x[k] <<"  "<< EE_y[k] <<"  "<< EE 
                  <<"  "<< Energy_Spectrum(double(k), Input_Parameters.i_Spectrum) << endl;
  }
  OutInitSpectrum.close();


  //dissi = 2.0*ens*vis;
  double u_p = sqrt(0.5*(int_x + int_y));
  double Taylor_scale = sqrt(TWO*u_p*u_p/(TWO*ens)) /*sqrt(15.0*u_p*u_p/(2.0*ens))*/;
  l_11 = PI*l_11/(2.0*u_p*u_p); 

  cout << "\n\n ==========================================================================\n"; 
  cout << " In spectral space:\n";
  cout << "\n <u^2> = " << int_x <<"  <v^2> = "<< int_y <<"  "
       << "u_rms = " << u_p 
       <<"\n L_11 = "<< l_11 << "  ens = "<< ens 
       << "\n Taylor_scale = " << Taylor_scale << endl;
  cout << " ==========================================================================" << endl;
  
  // Complex to real Fourier transforms
  physical = fftw_plan_dft_c2r_2d(Nx, Ny, uu, u,  FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);

  physical = fftw_plan_dft_c2r_2d(Nx, Ny, vv, v,  FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);

  
 
  // Output to a data file for initialization of the velocity flow field
  Write_Initial_Turbulence(QuadTree, Input_Parameters, u, v, Ny);

  //Write_Initial_Turbulence(InitMeshBlk, QuadTree, Input_Parameters, u, v, Ny);
  
  
  // Deallocations
  delete[] EE_x;
  delete[] EE_y;
  delete[] Crf;
  fftw_free(u);
  fftw_free(v);
  fftw_free(uu);
  fftw_free(vv);
   
}



// Compute the energy spectrum in Fourier space
void Compute_Spectrum(LESPremixed2D_Quad_Block *Soln_ptr,
                      AdaptiveBlock2D_List &Soln_Block_List,
                      LESPremixed2D_Input_Parameters &Input_Parameters,
                      double &u_average, double &v_average, double &vis) {

  double Lx, Ly, u_ave, v_ave;
  int Nx, Ny, ny;
 

  Lx = Input_Parameters.Box_Width;
  Ly = Input_Parameters.Box_Height;
  Nx = Input_Parameters.Number_of_Cells_Idir;
  Ny = Input_Parameters.Number_of_Cells_Jdir;
  ny = Ny/2+1;

  double        scale = 1.0/(Nx*Ny);  // Scale factor for the real to complex transform
  double        *u, *v;               // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv;             // Arrays to store the velocity fluctuations in Fourier space            
  fftw_plan     spectral;


  // Allocation of arrays
  u = (double *) malloc(Nx*Ny * sizeof(double));
  v = (double *) malloc(Nx*Ny * sizeof(double));
  uu = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
  vv = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
      
      
  // Mean velocities
  u_ave = u_average;
  v_ave = v_average;
    
  // Read and assign turbulent fluctuations from files for all the blocks
  Read_Turbulent_Fluctuations(Input_Parameters, u, v, u_ave, v_ave, Ny, scale);

  // Real to complex Fourier transform
  spectral = fftw_plan_dft_r2c_2d(Nx, Ny, u, uu, FFTW_ESTIMATE);
  fftw_execute(spectral); 
  fftw_destroy_plan(spectral);

  spectral = fftw_plan_dft_r2c_2d(Nx, Ny, v, vv,  FFTW_ESTIMATE);
  fftw_execute(spectral); 
  fftw_destroy_plan(spectral);


  double k1, k2;                    // wave numbers
  double *EE_x, *EE_y;              // spectra  
  ofstream OutSpectrum("Spectrum.dat"); // output file
  
      
  // Allocation of spectra
  double k_max = Max_k(Input_Parameters);
  cout <<"\n kmax = "<< k_max << endl;
  EE_x = new double[int(k_max)+1];
  EE_y = new double[int(k_max)+1];
  
  for (int p=0; p<=int(k_max); p++) {
    EE_x[p] = 0.0;
    EE_y[p] = 0.0;
  }
    
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      double abs_k = sqrt(k1*k1 + k2*k2);
      if (j == 0  ||  j == Ny/2) {
        EE_x[int(abs_k)] += 1.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 1.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      } else {
        EE_x[int(abs_k)] += 2.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 2.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      }
    }
  }
  
  double int_x = 0.0, int_y = 0.0; // spectra checks (integration of the spectra)
  double dissi, ens =0.0;          // dissipation of k in spectral space
  double l_11 = 0.0;               // L_11 integral length scale
  double EE, l;
 
  
  // Output the spectra to a file
  for (int k=0; k<=int(k_max); k++) {
    int_x += EE_x[k];
    int_y += EE_y[k];
    EE = 0.5*(EE_x[k]+EE_y[k]);
    if (k==0) {
      l_11 += 0.0; 
    } else {
      l_11 += EE/double(k);
    }  
    ens += EE*double(k*k);
    OutSpectrum << k <<"  "<< EE_x[k] <<"  "<< EE_y[k] <<"  "<< EE << "  "
            << double(k*k*k)*EE << endl;
  }
  
  dissi = 2.0*ens*vis;
  double u_p = sqrt(0.5*(int_x + int_y));
  double Taylor_scale = sqrt(TWO*u_p*u_p/(TWO*ens)) /*sqrt(15.0*u_p*u_p/(2.0*ens))*/;
  double Re_t = u_p*Taylor_scale/vis;
  double Kolmogorov_scale = pow(pow(vis, THREE)/dissi, 0.25);

  l_11 = PI*l_11/(2.0*u_p*u_p);
  l = 0.42*pow(u_p, 3.0)/dissi;
  
  cout << "\n\n ==========================================================================\n"; 
  cout << " In spectral space:\n";
  cout << "\n <u^2> = " << int_x <<"  <v^2> = "<< int_y <<"  "
       << "u_rms = " << u_p 
       << "\n L_11 = "<< l_11 <<"  eps = " << dissi << "  ens = "<< ens 
       << "\n l = " << l <<"  "<< "vis = " << vis << "  Re_t = " << Re_t
       << "\n Taylor scale = " << Taylor_scale 
       << "  Kolmogorov_scale = " << Kolmogorov_scale << endl;
  cout << " ==========================================================================" << endl;

  
  // Deallocations
  delete[] EE_x;
  delete[] EE_y;
  fftw_free(u);
  fftw_free(v);
  fftw_free(uu);
  fftw_free(vv);

}



void Rescale_Velocity(LESPremixed2D_Quad_Block *Soln_ptr,
                      AdaptiveBlock2D_List &Soln_Block_List,
                      LESPremixed2D_Input_Parameters &Input_Parameters,
                      const double &u_average, const double &v_average,
                      double &ko) {

      
  double k, Crf, u_ave, v_ave;
 
  k = Total_TKE(Soln_ptr, Soln_Block_List, Input_Parameters);
  Crf = sqrt(ko/k);
  //PRINT(ko); PRINT(k);  
  u_ave = u_average;
  v_ave = v_average;  

  for (int p = 0; p < Soln_Block_List.Nblk; ++p) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) { 
      for (int i = Soln_ptr[p].ICl; i <= Soln_ptr[p].ICu; ++i) {
        for (int j  = Soln_ptr[p].JCl; j <= Soln_ptr[p].JCu; ++j) {
          Soln_ptr[p].W[i][j].v.x = Crf*(Soln_ptr[p].W[i][j].v.x - u_ave); // + u_ave
          Soln_ptr[p].W[i][j].v.y = Crf*(Soln_ptr[p].W[i][j].v.y - v_ave); // + v_ave;
 	}
      }
    }
  }
  
  // Output the rescaled velocity field to a data file 
  Write_Rescaled_Velocity(Soln_ptr, Soln_Block_List, Input_Parameters);
  
}




void Rescale_Spectrum(// LESPremixed2D_Quad_Block *Soln_ptr,
//                       AdaptiveBlock2D_List &Soln_Block_List,
                      LESPremixed2D_Input_Parameters &Input_Parameters,
                      double &u_average, double &v_average, double &vis) {

  double Lx, Ly, u_ave, v_ave;
  int Nx, Ny, ny;
 

  Lx = Input_Parameters.Box_Width;
  Ly = Input_Parameters.Box_Height;
  Nx = Input_Parameters.Number_of_Cells_Idir;
  Ny = Input_Parameters.Number_of_Cells_Jdir;
  ny = Ny/2+1;

  double        scale = 1.0/(Nx*Ny);  // Scale factor for the real to complex transform
  double          *u, *v;             // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv;             // Arrays to store the velocity fluctuations in Fourier space
  fftw_plan    spectral, physical;
  

  // Allocation of arrays
  u = (double *) malloc(Nx*Ny * sizeof(double));
  v = (double *) malloc(Nx*Ny * sizeof(double));
  uu = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
  vv = (fftw_complex *) malloc(Nx*ny * sizeof(fftw_complex));
    
      
  // Mean velocities
  u_ave = u_average;
  v_ave = v_average;
    
  // Read and assign turbulent fluctuations from files for all the blocks
  Read_Turbulent_Fluctuations(Input_Parameters, u, v, u_ave, v_ave, Ny, scale);
  
  // Real to complex Fourier transforms
  spectral = fftw_plan_dft_r2c_2d(Nx, Ny, u, uu, FFTW_ESTIMATE);
  fftw_execute(spectral); 
  fftw_destroy_plan(spectral);

  spectral = fftw_plan_dft_r2c_2d(Nx, Ny, v, vv, FFTW_ESTIMATE);
  fftw_execute(spectral); 
  fftw_destroy_plan(spectral);

    
  
  double k1, k2, abs_k;                          // wave numbers
  double *EE_x, *EE_y, *Crf;                     // spectra and rescaling factor  
  ofstream OutNewSpectrum("Rescaled_Spectrum.dat"); // output file
  
      
  // Allocation of spectra
  double k_max = Max_k(Input_Parameters);
  cout <<"\n kmax = "<< k_max << endl;
  EE_x = new double[int(k_max)+1];
  EE_y = new double[int(k_max)+1];
  Crf  = new double[int(k_max)+1];
  
  for (int p=0; p<=int(k_max); p++) {
    EE_x[p] = 0.0;
    EE_y[p] = 0.0;
  }
  

  //------------------  Calculate the actual spectrum ---------------------//  
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      if (j == 0  ||  j == Ny/2) {
        EE_x[int(abs_k)] += 1.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 1.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      } else {
        EE_x[int(abs_k)] += 2.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 2.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      }
    }
  }
  
  double int_x = 0.0, int_y = 0.0; // spectra checks (integration of the spectra)
  double dissi, ens = 0.0;         // dissipation of k in spectral space
  double l_11 = 0.0;               // L_11 integral length scale
  double EE, Eo, l;
 
  
  //------ Calculate the scaling factor of the velocity fluctuations  -------//
  for (int k=0; k<=int(k_max); k++) {
    // int_x += EE_x[k];
//     int_y += EE_y[k];
    EE = 0.5*(EE_x[k]+EE_y[k]);
    Eo = Energy_Spectrum(double(k), Input_Parameters.i_Spectrum);
    Crf[k] = sqrt(Eo/EE);
  }
  

  //-----------      Rescale fluctuations       -------------//
  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      // Rescaling
      uu[i*ny+j][0] = Crf[int(abs_k)]*uu[i*ny+j][0];
      uu[i*ny+j][1] = Crf[int(abs_k)]*uu[i*ny+j][1]; 
      vv[i*ny+j][0] = Crf[int(abs_k)]*vv[i*ny+j][0];
      vv[i*ny+j][1] = Crf[int(abs_k)]*vv[i*ny+j][1];   
    }
  }


  //----------    Calculate the new (rescaled) energy spectrum   ---------//

  // Reset EE_x and EE_y to ZERO
  for (int p=0; p<=int(k_max); p++) {
    EE_x[p] = 0.0;
    EE_y[p] = 0.0;
  }

  for (int i=0; i<Nx; i++) {
    for (int j=0; j<Ny/2+1; j++) {
      // Components of the wave number vector
      if (i<=Nx/2) {
        k1 = k_1(i);
      } else {
        k1 = k_1(i-Nx);
      }
      k2 = k_2(j);
      // Wave number magnitude
      abs_k = sqrt(k1*k1 + k2*k2);
      if (j == 0  ||  j == Ny/2) {
        EE_x[int(abs_k)] += 1.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 1.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      } else {
        EE_x[int(abs_k)] += 2.0*( sqr(uu[i*ny+j][0]) + sqr(uu[i*ny+j][1]) );
        EE_y[int(abs_k)] += 2.0*( sqr(vv[i*ny+j][0]) + sqr(vv[i*ny+j][1]) );
      }
    }
  }

 
  //------------    Output the new spectrum to a file   --------------//
  int_x = ZERO;  
  int_y = ZERO;
  for (int k=0; k<=int(k_max); k++) {
    int_x += EE_x[k];
    int_y += EE_y[k];
    EE = 0.5*(EE_x[k]+EE_y[k]);
    if (k==0) {
      l_11 += 0.0; 
    } else {
      l_11 += EE/double(k);
    }  
    ens += EE*double(k*k);
    OutNewSpectrum << k <<"  "<< EE_x[k] <<"  "<< EE_y[k] <<"  "<< EE << "  "
                   << double(k*k*k)*EE << "  " <<  Energy_Spectrum(double(k), Input_Parameters.i_Spectrum)
                   << "  " << Crf[k] << endl;
  }

  dissi = 2.0*ens*vis;
  double u_p = sqrt(0.5*(int_x + int_y));
  double Taylor_scale = sqrt(TWO*u_p*u_p/(TWO*ens)) /*sqrt(15.0*u_p*u_p/(2.0*ens))*/;
  double Re_t = u_p*Taylor_scale/vis;
  double Kolmogorov_scale = pow(pow(vis, THREE)/dissi, 0.25);
 
  l_11 = PI*l_11/(2.0*u_p*u_p);
  l = 0.42*pow(u_p, 3.0)/dissi;
  
  cout << "\n\n ==========================================================================\n"; 
  cout << " In spectral space for the rescaled field:\n";
  cout << "\n <u^2> = " << int_x <<"  <v^2> = "<< int_y <<"  "
       << "u_rms = " << u_p 
       << "\n L_11 = "<< l_11 <<"  eps = " << dissi << "  ens = "<< ens 
       << "\n l = " << l <<"  "<< "vis = " << vis << "  Re_t = " << Re_t
       << "\n Taylor scale = " << Taylor_scale 
       << "  Kolmogorov_scale = " << Kolmogorov_scale << endl;
  cout << " ==========================================================================" << endl;
  //  cout << "\n In spectral space for the rescaled field:" << endl;
//   cout << " <u^2> = " << int_x <<"  <v^2> = "<< int_y <<"  "
//        << "u_rms = " << u_p <<"\n L_11 = "<< l_11 <<"  eps = " << dissi 
//        << "  ens = "<< ens << "  l = " << l <<"  "<< "vis = " << vis 
//        << "\n Taylor scale = " << Taylor_scale 
//        << "  Kolmogorov_scale = " << Kolmogorov_scale 
//        << "  Re_t = " << Re_t << endl;


  // Complex to real Fourier transforms
  physical = fftw_plan_dft_c2r_2d(Nx, Ny, uu, u,  FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);

  physical = fftw_plan_dft_c2r_2d(Nx, Ny, vv, v,  FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);


  // Output the rescaled velocity field to a data file 
  Write_Rescaled_Velocity(Input_Parameters, u, v, Ny);

  
  // Deallocations
  delete[] EE_x;
  delete[] EE_y;
  delete[] Crf;
  fftw_free(u);
  fftw_free(v);
  fftw_free(uu);
  fftw_free(vv);

}







/*================================================*\
            Write/Read file tools        
\*================================================*/

// Write the initial turbulent velocity fluctuations for the initialization
// of a homogenous isotropic turbulent flow field
void Write_Initial_Turbulence(QuadTreeBlock_DataStructure &QuadTree,
			      LESPremixed2D_Input_Parameters &Input_Parameters,
                              double *u, double *v, const int &Ny) {

  int Nblks_i, Nblks_j, ICblk, JCblk, ii, jj;
  char blk_num[10];
  string block, filename;
  int Block_Number;
  ofstream turbulencedata_out;

  Nblks_i = Input_Parameters.Number_of_Blocks_Idir;
  Nblks_j = Input_Parameters.Number_of_Blocks_Jdir;  
  ICblk = Input_Parameters.Number_of_Cells_Idir/Nblks_i;  // Number of cells per block in the I direction
  JCblk = Input_Parameters.Number_of_Cells_Jdir/Nblks_j;  // Number of cells per block in the J direction
  

  // turbulencedata_out.open("Initial_Turbulence.dat", ios::out);
//   if(turbulencedata_out.fail()){
//     cerr<<"\nError opening file: Initial_Turbulence.dat to write" << endl;
//     exit(1);
//   }
  
  
  for (int p=0; p<Nblks_i; ++p) {
    for (int q=0; q<Nblks_j; ++q) {
      filename = "Initial_Turbulence_";
      block = "Block_"; 
      Block_Number = QuadTree.Roots[p][q].block.gblknum; //q*Nblks_i + p;
      sprintf(blk_num, "%.6d", Block_Number);
      block += blk_num;
      filename += block + ".dat";

      // Open file
      turbulencedata_out.open(filename.c_str(), ios::out);
      if (turbulencedata_out.fail()) {
	cerr<<"\nError opening file: " << filename <<" to write." << endl;
	exit(1);
      }

      //turbulencedata_out << block << endl;
      turbulencedata_out.setf(ios::scientific); 
      for (int i=0; i<ICblk; ++i) {
        ii = p*ICblk + i;
        for (int j=0; j<JCblk; ++j) {
          jj = q*JCblk + j;
          turbulencedata_out << setprecision(10) << u[ii*Ny+jj] << " " << v[ii*Ny+jj] << "\n";
        }
      }
      turbulencedata_out.unsetf(ios::scientific);
      turbulencedata_out.close();

      //turbulencedata_out << endl;
    }
  }

  //turbulencedata_out.close();   

}


////////////////////  OVERLOADED  //////////////////////////////////////////
// Write the initial turbulent velocity fluctuations for the initialization
// of a homogenous isotropic turbulent flow field
void Write_Initial_Turbulence(Grid2D_Quad_Block   **InitMeshBlk,
			      QuadTreeBlock_DataStructure &QuadTree,
			      LESPremixed2D_Input_Parameters &Input_Parameters,
                              double *u, double *v, 
			      const int &Ny) {

  int Nblks_i, Nblks_j, ICblk, JCblk, ii, jj;
  char blk_num[10];
  string block, filename;
  int Block_Number;
  ofstream turbulencedata_out;

  Nblks_i = Input_Parameters.Number_of_Blocks_Idir;
  Nblks_j = Input_Parameters.Number_of_Blocks_Jdir;  
  ICblk = Input_Parameters.Number_of_Cells_Idir/Nblks_i;  // Number of cells per block in the I direction
  JCblk = Input_Parameters.Number_of_Cells_Jdir/Nblks_j;  // Number of cells per block in the J direction
  

  turbulencedata_out.open("Initial_Turbulence.dat", ios::out);
  if(turbulencedata_out.fail()){
    cerr<<"\nError opening file: Initial_Turbulence.dat to write" << endl;
    exit(1);
  }
  
  
  for (int i_blk=0; i_blk<Nblks_i; ++i_blk) {
    for (int j_blk=0; j_blk<Nblks_j; ++j_blk) {
      if ( InitMeshBlk[i_blk][j_blk].Node != NULL ) { // Mesh block is used!!!!

	int ICl = InitMeshBlk[i_blk][j_blk].ICl, JCl = InitMeshBlk[i_blk][j_blk].JCl; 
	int ICu = InitMeshBlk[i_blk][j_blk].ICu, JCu = InitMeshBlk[i_blk][j_blk].JCu;

	if (ICu-ICl+1 != ICblk  ||  JCu-JCl+1 != JCblk) {
	  cout << "\nERROR: Number of cell does not match writing Initial_Turbulence.dat";
	}
		
	// PRINT(QuadTree.Roots[i_blk][j_blk].block.gblknum);
// 	PRINT(j_blk*Nblks_i + i_blk);

	turbulencedata_out.setf(ios::scientific); 
	for (int i=0; i<ICblk; ++i) {
	  ii = i_blk*ICblk + i;
	  for (int j=0; j<JCblk; ++j) {
	    jj = j_blk*JCblk + j;
	    turbulencedata_out << setprecision(10) <<  InitMeshBlk[i_blk][j_blk].Cell[i+ICl][j+JCl].Xc  
			       <<  " " << u[ii*Ny+jj] << " " << v[ii*Ny+jj] << "\n";
	  }
	}
	turbulencedata_out.unsetf(ios::scientific);
	turbulencedata_out << endl;

      } // end if
    }
  }
  turbulencedata_out.close();   

}


 


void Write_Rescaled_Velocity(LESPremixed2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LESPremixed2D_Input_Parameters &Input_Parameters) {

  char blk_num[10];
  string block, filename;
  ofstream out_file;


  /* Write the rescaled velocity for each solution block. */ 
  for (int p = 0 ; p <= Soln_Block_List.Nblk-1 ; ++p ) {
    if (Soln_Block_List.Block[p].used == ADAPTIVEBLOCK2D_USED) {
      // File name based on global block number.
      filename = "Rescaled_Velocity_";
      block = "Block_";
      sprintf(blk_num, "%.6d", Soln_Block_List.Block[p].gblknum);
      block += blk_num;
      filename += block + ".dat";
    
      // Open file
      out_file.open(filename.c_str(), ios::out);
      if (out_file.fail()) {
	cerr<<"\nError opening file: " << filename <<" to write." << endl;
	exit(1);
      }

      //out_file << block << endl;
      out_file.setf(ios::scientific);
      for (int i = Soln_ptr[p].ICl/*-Soln_ptr[p].Nghost*/; i <= Soln_ptr[p].ICu/*+Soln_ptr[p].Nghost*/; ++i) {
	for (int j = Soln_ptr[p].JCl/*-Soln_ptr[p].Nghost*/; j <= Soln_ptr[p].JCu/*+Soln_ptr[p].Nghost*/; ++j) {
	  out_file << setprecision(14) << Soln_ptr[p].W[i][j].v.x <<" "<< Soln_ptr[p].W[i][j].v.y << endl;
	}
      }
      out_file.unsetf(ios::scientific);
      // Close file.
      out_file.close(); 

    }
  }

}



void Write_Rescaled_Velocity(LESPremixed2D_Input_Parameters &Input_Parameters,
                             double *u, double *v, const int &Ny) {

  int Nblks_i, Nblks_j, ICblk, JCblk, ii, jj;
  char blk_num[10];
  string block;
  int Block_Number;
  ofstream out_rescaled;

  Nblks_i = Input_Parameters.Number_of_Blocks_Idir;
  Nblks_j = Input_Parameters.Number_of_Blocks_Jdir;  
  ICblk = Input_Parameters.Number_of_Cells_Idir/Nblks_i;  // Number of cells per block in the I direction
  JCblk = Input_Parameters.Number_of_Cells_Jdir/Nblks_j;  // Number of cells per block in the J direction
  

  out_rescaled.open("Rescaled_Velocity.dat", ios::out);
  if(out_rescaled.fail()){
    cerr<<"\nError opening file: Rescaled_Velocity.dat to write" << endl;
    exit(1);
  }
  
  
  for (int p=0; p<Nblks_i; ++p) {
    for (int q=0; q<Nblks_j; ++q) {
      block = "Block_"; 
      Block_Number = q*Nblks_i + p;
      sprintf(blk_num, "%.6d", Block_Number);
      block += blk_num;
      out_rescaled << block << endl;
      out_rescaled.setf(ios::scientific); 
      for (int i=0; i<ICblk; ++i) {
        ii = p*ICblk + i;
        for (int j=0; j<JCblk; ++j) {
          jj = q*JCblk + j;
          out_rescaled << setprecision(10) << u[ii*Ny+jj] << " " << v[ii*Ny+jj] << "\n";
        }
      }
      out_rescaled.unsetf(ios::scientific);
      out_rescaled << endl;
    }
  }
  out_rescaled.close();   


}



// Read the turbulent velocity fluctuations for computing the energy spectrum
void Read_Turbulent_Fluctuations(LESPremixed2D_Input_Parameters &Input_Parameters,
                                 double *u, double *v, 
                                 double &u_ave, double &v_ave,
                                 int &Ny, double &scale) {

  string block, file;
  char blk_num[10];
  ifstream turbulence_in;
  int ii, jj, Nblks_i, Nblks_j, ICblk, JCblk, Blknum;

  Nblks_i = Input_Parameters.Number_of_Blocks_Idir;
  Nblks_j = Input_Parameters.Number_of_Blocks_Jdir;  
  ICblk = Input_Parameters.Number_of_Cells_Idir/Nblks_i;
  JCblk = Input_Parameters.Number_of_Cells_Jdir/Nblks_j;

  for (int p=0; p<Nblks_i; p++) {
    for (int q=0; q<Nblks_j; q++) {
      Blknum = q*Nblks_i+p;
      file = "Initial_Turbulence_";
      block = "Blk_";
      sprintf(blk_num, "%.6d", Blknum);
      block += blk_num;
      file += block + ".dat";
  
      // Open turbulence data file for reading
      turbulence_in.open(file.c_str(), ios::in); 
      // Check to see if successful
      if (turbulence_in.fail()){ 
        cerr<<"\nError opening file: "<< file <<" to read" <<endl;
        exit(1); 
      } 
              
      turbulence_in.setf(ios::skipws);     
      for (int i = 0; i < ICblk; i++ ) {
        ii = p*ICblk + i;
        for (int j = 0; j < JCblk; j++ ) {
          jj = q*JCblk + j;
          turbulence_in >> u[ii*Ny+jj] >> v[ii*Ny+jj];
          u[ii*Ny+jj] = scale*(u[ii*Ny+jj] - u_ave);
          v[ii*Ny+jj] = scale*(v[ii*Ny+jj] - v_ave);          
	}
      }
      turbulence_in.unsetf(ios::skipws);
      turbulence_in.close();
      //cout <<"\nDone with " << block <<endl;    
    }
  }

}


/*-------------------------------------------------------*\
  Calculate the longitudinal integral length scale L11
\*-------------------------------------------------------*/
int Longitudinal_Correlation(QuadTreeBlock_DataStructure &QuadTree,
                             AdaptiveBlockResourceList  &Global_Soln_Block_List,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LESPremixed2D_Quad_Block  *Soln_ptr,   
			     LESPremixed2D_Input_Parameters  &Input_Parameters, 
			     const double & u_ave, const double &sqr_u) {

  
  ofstream Out_Corr_Function;
  int error_flag = 0;
  int Nblks = Global_Soln_Block_List.Nused;
  double xref, yref, uref, ucorr;
  double Yfuel_conditional = ZERO;

  
  
 
  int gblknum;  
  double r, dr, max_r, fr, Lx, area, local_A, total_A, R11, L11;
  Lx = Input_Parameters.Box_Width; 
  dr = MILLION;

  //Conditional longitudinal correlation on fresh gas
  if (Input_Parameters.react_name != "NO_REACTIONS") {
    Yfuel_conditional = 0.95*Input_Parameters.Fresh_Fuel_Mass_Fraction;
  }

  for (int k = 0; k < Soln_Block_List.Nblk; ++k) {
    if (Soln_Block_List.Block[k].used == ADAPTIVEBLOCK2D_USED) {
      for (int ii = Soln_ptr[k].ICl; ii <= Soln_ptr[k].ICu; ++ii) {
        for (int jj = Soln_ptr[k].JCl; jj <= Soln_ptr[k].JCu; ++jj) {
	  if (Soln_ptr[k].W[ii][jj].spec[0].c >= Yfuel_conditional) {
	    dr = min(Soln_ptr[k].Grid.Cell[ii][jj].Xc.x - Soln_ptr[k].Grid.Cell[ii-1][jj].Xc.x, dr);          
	  }
	}
      }
    }
  }
 
  dr = CFFC_Minimum_MPI(dr);  

  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.open("Correlation_Function.dat", ios::out);
    if(Out_Corr_Function.fail()){
      cerr<<"\nError opening file: Correlation_Function.dat to write" << endl;
      exit(1);
    }
  }
  
  int *new_blocks_CPU, *CPUs_in_new_blocks;
  int my_rank, undefined_rank, new_blocks_BLK;
  CPUs_in_new_blocks = new int[Global_Soln_Block_List.Ncpu];
  new_blocks_CPU = new int[Global_Soln_Block_List.Ncpu];
  
#ifdef _MPI_VERSION
    MPI::Intracomm new_comm;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     new_group;
    undefined_rank = MPI::UNDEFINED;
#else
    undefined_rank = -1;
#endif
  
  LESPremixed2D_Quad_Block  SolnBlk_duplicated;
  QuadTreeBlock *quadtree_block_duplicated_ptr;
  int SolnBlk_duplicated_info_level;

  Vector2D Vcorr, Xcorr;
  Vcorr.zero();  Xcorr.zero();
  bool correlated_flag, flag = false;
  int count = 0, count1 = 0;
  r = ZERO;
  L11 = ZERO;

  if (Input_Parameters.i_Grid == GRID_PERIODIC_BOX &&
      Input_Parameters.react_name == "NO_REACTIONS") {
    max_r = HALF*(Lx-dr);
  } else  {
    max_r = Lx-dr;
  }
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  while( r <= max_r) {
    total_A = ZERO;
    R11 = ZERO;
    correlated_flag = false;
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.  

    for (int iCPU = 0; iCPU <= QuadTree.Ncpu-1; ++iCPU ) { // Loop over available processors.
      for (int iBLK = 0; iBLK <= QuadTree.Nblk-1; ++iBLK ) { // Loop over available blocks.
        if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
          if (QuadTree.Blocks[iCPU][iBLK]->block.used) { // Check if the solution block is used.
        
            quadtree_block_duplicated_ptr = QuadTree.Blocks[iCPU][iBLK];
            new_blocks_CPU[0] = quadtree_block_duplicated_ptr->block.info.cpu;
            new_blocks_BLK = quadtree_block_duplicated_ptr->block.info.blknum;
            CPUs_in_new_blocks[0] = new_blocks_CPU[0];

	    if (Global_Soln_Block_List.Ncpu > 1) {
	      for (int iNEW = 1 ; iNEW < Global_Soln_Block_List.Ncpu; ++iNEW ) {
		if (iNEW != new_blocks_CPU[0]) {
		  new_blocks_CPU[iNEW] = iNEW;
		} else {
		  new_blocks_CPU[iNEW] = 0;
		}  
		CPUs_in_new_blocks[iNEW] = new_blocks_CPU[iNEW];
	      }
	    }
  
#ifdef _MPI_VERSION
	    new_group = big_group.Incl(Global_Soln_Block_List.Ncpu, CPUs_in_new_blocks);
	    new_comm = MPI::COMM_WORLD.Create(new_group);
#endif
	    if (Soln_Block_List.ThisCPU == new_blocks_CPU[0]) {  
	      Copy_Solution_Block(SolnBlk_duplicated, Soln_ptr[new_blocks_BLK]);
	      SolnBlk_duplicated_info_level = Soln_Block_List.Block[new_blocks_BLK].info.level;
	    }	   
#ifdef _MPI_VERSION
	    if (my_rank != undefined_rank) { 
	      Broadcast_Solution_Block(SolnBlk_duplicated, new_comm, new_blocks_CPU[0]);
	      new_comm.Bcast(&SolnBlk_duplicated_info_level, 1, MPI::INT, 0);
	    }
#endif

#ifdef _MPI_VERSION
	    if (new_comm != MPI::COMM_NULL) new_comm.Free();
	    new_group.Free();
#endif
	    for (int i_ref = SolnBlk_duplicated.ICl; i_ref <= SolnBlk_duplicated.ICu; ++i_ref) {
	      for (int j_ref = SolnBlk_duplicated.JCl; j_ref <= SolnBlk_duplicated.JCu; ++j_ref) {
                if (SolnBlk_duplicated.W[i_ref][j_ref].spec[0].c >= Yfuel_conditional) {
		  xref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref].Xc.x;
		  yref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref].Xc.y;
		  uref = SolnBlk_duplicated.W[i_ref][j_ref].v.x;
		  flag = false;
		} else {
                  continue;
		}
	       
		for (int q = 0; q < Soln_Block_List.Nblk; ++q) {
		  if (Soln_Block_List.Block[q].used == ADAPTIVEBLOCK2D_USED) {
		    for (int i = Soln_ptr[q].ICl; i <= Soln_ptr[q].ICu; ++i) {
		      for (int j = Soln_ptr[q].JCl; j <= Soln_ptr[q].JCu; ++j) {

			if ((Soln_ptr[q].W[i][j].spec[0].c >= Yfuel_conditional) &&
                            (Soln_ptr[q].Grid.xfaceW(i,j).x < (xref + r) &&
			     Soln_ptr[q].Grid.xfaceE(i,j).x > (xref + r))) {

			  // Finer reference block than local block
			  if (SolnBlk_duplicated_info_level >= Soln_Block_List.Block[q].info.level  &&
			      (Soln_ptr[q].Grid.xfaceS(i,j).y < yref  &&  Soln_ptr[q].Grid.xfaceN(i,j).y > yref)) {

			    if (Soln_ptr[q].Grid.Cell[i][j].Xc.y == yref &&
				Soln_ptr[q].Grid.Cell[i][j].Xc.x == (xref + r)) {
			      Vcorr.x = Soln_ptr[q].W[i][j].v.x;
			    } else {
			      Vector2D dX;
			      dX.x = (xref + r) - Soln_ptr[q].Grid.Cell[i][j].Xc.x;
			      dX.y = yref - Soln_ptr[q].Grid.Cell[i][j].Xc.y;
			      Vcorr.x = Soln_ptr[q].W[i][j].v.x + Soln_ptr[q].phi[i][j].v.x 
				* (Soln_ptr[q].dWdx[i][j].v.x * dX.x + Soln_ptr[q].dWdy[i][j].v.x * dX.y);
			    }
			    // correlated u
			    ucorr = Vcorr.x;
			    local_A = Soln_ptr[q].Grid.Cell[i][j].A;
			    count1 ++;
			    flag = true;

			  // Coarser reference block than local block  
			  } else if (SolnBlk_duplicated_info_level < Soln_Block_List.Block[q].info.level  &&
				     (Soln_ptr[q].Grid.Cell[i][j].Xc.y < yref  &&  Soln_ptr[q].Grid.Cell[i][j+1].Xc.y > yref)) {
			    Xcorr.x = xref + r;
			    Xcorr.y = yref; 
			    Bilinear_Interpolation(Soln_ptr[q].W[i][j].v, Soln_ptr[q].Grid.Cell[i][j].Xc,
						   Soln_ptr[q].W[i][j+1].v, Soln_ptr[q].Grid.Cell[i][j+1].Xc,
						   Soln_ptr[q].W[i+1][j+1].v, Soln_ptr[q].Grid.Cell[i+1][j+1].Xc,
						   Soln_ptr[q].W[i+1][j].v, Soln_ptr[q].Grid.Cell[i+1][j].Xc,
						   Xcorr, Vcorr);
			    // correlated u
			    ucorr = Vcorr.x;
			    local_A = Soln_ptr[q].Grid.Cell[i][j].A;
			    count1 ++;
			    flag = true;
			  } // end if
			
                                                                                               
// 			  if (Soln_ptr[q].Grid.Cell[i][j].Xc.y == yref &&
// 			      Soln_ptr[q].Grid.Cell[i][j].Xc.x == (xref + r)) {
// 			    Vcorr.x = Soln_ptr[q].W[i][j].v.x;
// 			  } else {
// 			    Xcorr.x = xref + r;
// 			    Xcorr.y = yref; 
// 			    // Bilinear_Interpolation(Soln_ptr[q].WnNW(i,j).v, Soln_ptr[q].Grid.nodeNW(i,j).X,
// 			    // 						 Soln_ptr[q].WnNE(i,j).v, Soln_ptr[q].Grid.nodeNE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSE(i,j).v, Soln_ptr[q].Grid.nodeSE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSW(i,j).v, Soln_ptr[q].Grid.nodeSW(i,j).X,
// 			    // 						 Xcorr, Vcorr);
// 			    Bilinear_Interpolation(Soln_ptr[q].W[i-1][j].v, Soln_ptr[q].Grid.Cell[i-1][j].Xc,
// 						   Soln_ptr[q].W[i][j+1].v, Soln_ptr[q].Grid.Cell[i][j+1].Xc,
// 						   Soln_ptr[q].W[i+1][j].v, Soln_ptr[q].Grid.Cell[i+1][j].Xc,
// 						   Soln_ptr[q].W[i][j-1].v, Soln_ptr[q].Grid.Cell[i][j-1].Xc,
// 						   Xcorr, Vcorr);
// 			  }

			  
			} // end if

		      }                       
		    }
		  } //end if
		  if (flag) break;
		} // end for
        
		if (flag) {
		  area = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref].A + local_A;
		  total_A += area;
		  R11 += (ucorr - u_ave)*(uref - u_ave)*area;
		  correlated_flag = true;
		  count++;
		}    
           
	      } // end for
	    } // end for
	 
	  }
	}
      }
    }
  
    if (CFFC_OR_MPI(correlated_flag)) {
      total_A = CFFC_Summation_MPI(total_A);
      R11 = CFFC_Summation_MPI(R11);
       
      R11 = R11/total_A;     // Area-averaged two-point correlation
      fr = R11/sqr_u;        // Longitudinal autocorrelation function
      L11 += fr*dr;          // Integrate the above funtion to obtain L11
      if (CFFC_Primary_MPI_Processor()) {
        Out_Corr_Function << r << "  " << fr << endl;  
      }
    }
       
    //cout << "\n->" << r; 
    r += dr;
  } // end while
  
  count = CFFC_Summation_MPI(count);
  count1 = CFFC_Summation_MPI(count1);
    
  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.close();
    cout << "\n\n *** L11 = " << L11 << " ***" << endl;
  }

  delete[] CPUs_in_new_blocks;   CPUs_in_new_blocks = NULL;
  delete[] new_blocks_CPU;       new_blocks_CPU = NULL;
  
  if (CPUs_in_new_blocks != NULL || new_blocks_CPU != NULL) error_flag = 1;
       
  return (error_flag);
}



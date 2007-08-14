#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif //MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif //CFD_INCLUDED

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "TurbulenceModelling.h"
#endif // TURBULENCE_MODELLING_INCLUDED   

// initalize the static members of the primitive and conservative classes.
double Turbulence_Model_k_omega::alpha = 0.52;
//13.0/25.0; // 5.0/9.0 (1989 model coefficient)
double Turbulence_Model_k_omega::sigma = 0.5;
double Turbulence_Model_k_omega::sigma_star = 0.5;
double Turbulence_Model_k_omega::beta = 0.072;
//9.0/125.0; //3.0/40.0 (1989 model coefficient)
double Turbulence_Model_k_omega::f_beta_const = 1.0;
double Turbulence_Model_k_omega::beta_star = 9.0/100.0;
double Turbulence_Model_k_omega::f_beta_star_const = 1.0;
double Turbulence_Model_k_omega::EDM_const = 1.0;
double Turbulence_Model_k_omega::y_sublayer = 2.5;
double Turbulence_Model_k_omega::Karman_const = 0.41;
double Turbulence_Model_k_omega::C_const = 5.0;
double Turbulence_Model_k_omega::yplus_l = 30.0; 
double Turbulence_Model_k_omega::yplus_u =250.0; 



// f_beta and f_betastar (1998 k-omega model coefficient)
// (Wilcox textbook P121)
double Turbulence_Model_k_omega::f_beta(const Tensor3D &rotation, 
                                        const Tensor3D &strain_rate,
                                        const double &omega)const{
   
   double kai_w, kai_w_numerator, base; 

//      |Omega_ij Omega_jk S_ki|
// Xw = |----------------------| 
//      |  (beta_o* omega)^3   |
   kai_w_numerator = rotation.xy*(-rotation.xy)*strain_rate.xx + 
      rotation.xy*rotation.yz*strain_rate.xz + 
      rotation.xz*(-rotation.xz)*strain_rate.xx +  
      rotation.xz*(-rotation.yz)*strain_rate.xy + 
      (-rotation.xy)*rotation.xy*strain_rate.yy + 
      (-rotation.xy)*rotation.xz*strain_rate.yz + 
      rotation.yz*(-rotation.xz)*strain_rate.xy +  
      rotation.yz*(-rotation.yz)*strain_rate.yy +
      (-rotation.xz)*rotation.xy*strain_rate.yz + 
      (-rotation.xz)*rotation.xz*strain_rate.zz + 
      (-rotation.yz)*(-rotation.xy)*strain_rate.xz +  
      (rotation.yz)*rotation.yz*strain_rate.zz ;
   
   base = beta_star*max(TOLER, omega);
   kai_w =fabs(kai_w_numerator/max(TOLER, pow(base, 3)));
   // cout<<" \n kw = "<<fabs((ONE+70*kai_w)/(ONE+80*kai_w));
   return ((ONE+70*kai_w)/(ONE+80*kai_w));
   

}
double Turbulence_Model_k_omega::f_betastar(
   const double &dkdx,const double &dkdy, const double &dkdz, 
   const double &domegadx, const double &domegady, 
   const double &domegadz, const double &omega)const{
   
   
   double kai_k,base;
//          gradient(k) gradient (omega)
//  Xk = ---------------------------------
//          omega^3
   base = max(TOLER, omega);
   kai_k = (dkdx*domegadx + dkdy*domegady + dkdz*domegadz)/max(TOLER, pow(base, 3));
   
   if (kai_k <= ZERO){
      //    cout<<"\n kai_k= "<<  kai_k <<"      ---      "<<ONE;
      return (1.0);
   }else{
      //  cout<<"\n kai_k= "<<  kai_k <<"      ---      "
      return ((1.0+680*sqr(kai_k))/(1.0+400*sqr(kai_k)));
   }
   
}

double Turbulence_Model_k_omega::omega_sublayer_BC(
   const double &d, const double &mu, const double &y) const {
   return (SIX*mu/(d*beta*y*y));
}

// wall function
int Turbulence_Model_k_omega::wall_function(
   const double &d, const double &mu, Turbulent3DWallData &wall , 
   double &dk, double &domega){

   double ke, omegae;
   
   if(wall.yplus <= yplus_l) {
     
      
      //  cout<<wall.yplus<<endl;
      
   //   wall.utau = 0.3;// xinfeng notes: this is theexperimental value
      
      // Set k
      ke = sqr(wall.utau)/sqrt(beta_star);
      // Set Omega
      if(wall.utau == ZERO){
         omegae = 6.0*(mu/d)/(beta*sqr(wall.ywall));
      }else{   
         omegae = wall.utau/(sqrt(beta_star)*Karman_const*wall.ywall);   
         
      } /* endif */
      
      // Update conserved values
      dk =  d*ke;
      domega = d*omegae;

       
      // } /* endif */
      
   } /* endif */
   
   return 0;
   
} //end of wall function

int Turbulence_Model_k_omega::low_Reynolds_number_formulation(
   const double &d, const double &mu, 
   const int i, const int j, const int k,
   Turbulent3DWallData &wall,
   Grid3D_Hexa_Block &grid,
   double &domega){
   
   
   if ((wall.yplus <= y_sublayer)
       ||((i== grid.ICl) &&(grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeW[j][k] == BC_NO_SLIP  ||
                             grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeW[j][k] == BC_ADIABATIC_WALL)) 
       ||((i==grid.ICu) &&(grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
                            grid.BCtypeE[j][k] == BC_NO_SLIP  ||
                            grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeE[j][k] == BC_ADIABATIC_WALL))
       ||((j == grid.JCl) &&(grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
                              grid.BCtypeS[i][k] == BC_NO_SLIP  ||
                              grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                              grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                              grid.BCtypeS[i][k] == BC_ADIABATIC_WALL))
       ||((j ==grid.JCu) &&(grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeN[i][k] == BC_NO_SLIP  ||
                             grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeN[i][k] == BC_ADIABATIC_WALL))  
       ||((k == grid.KCl) &&(grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
                              grid.BCtypeB[i][j] == BC_NO_SLIP  ||
                              grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                              grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                              grid.BCtypeB[i][j] == BC_ADIABATIC_WALL)) 
       ||((k ==grid.KCu) &&(grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
                             grid.BCtypeT[i][j] == BC_NO_SLIP  ||
                             grid.BCtypeT[i][j]  == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeT[i][j] == BC_ADIABATIC_WALL))){
      
      domega =  d*6.0*(mu/d)/(beta*sqr(wall.ywall));
      
        
   }// end of checking
    
   return 0;
      
   
}
int Turbulence_Model_k_omega::wall_function_residual_evaluation(
   Turbulent3DWallData &wall, double &dUdt_rhok, double &dUdt_rhoomega){
   
   // don't need to update the log-layer soulution for rhok and rhoomega  
   if(wall.yplus <= yplus_u) {      
      dUdt_rhok =  0.0;
      dUdt_rhoomega = 0.0;
      
   } /* endif */
   
   return 0;
     
}

int Turbulence_Model_k_omega::low_Reynolds_number_formulation_residual_evaluaton(
   const int i, const int j, const int k, 
   Turbulent3DWallData &wall,   Grid3D_Hexa_Block &grid,
   double &dUdt_rhoomega){

   if ((wall.yplus <= y_sublayer)
       ||((i== grid.ICl) &&(grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeW[j][k] == BC_NO_SLIP  ||
                             grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeW[j][k] == BC_ADIABATIC_WALL)) 
       ||((i==grid.ICu) &&(grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
                            grid.BCtypeE[j][k] == BC_NO_SLIP  ||
                            grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeE[j][k] == BC_ADIABATIC_WALL))
       ||((j == grid.JCl) &&(grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
                              grid.BCtypeS[i][k] == BC_NO_SLIP  ||
                              grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                              grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                              grid.BCtypeS[i][k] == BC_ADIABATIC_WALL))
       ||((j ==grid.JCu) &&(grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeN[i][k] == BC_NO_SLIP  ||
                             grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeN[i][k] == BC_ADIABATIC_WALL))  
       ||((k == grid.KCl) &&(grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
                              grid.BCtypeB[i][j] == BC_NO_SLIP  ||
                              grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                              grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                              grid.BCtypeB[i][j] == BC_ADIABATIC_WALL)) 
       ||((k ==grid.KCu) &&(grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
                             grid.BCtypeT[i][j] == BC_NO_SLIP  ||
                             grid.BCtypeT[i][j]  == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeT[i][j] == BC_ADIABATIC_WALL))){
      
      dUdt_rhoomega =  0.0;
    
   }// end of checking
    
   return 0;
      


   
}


int Turbulence_Model_k_omega::automatic_wall_treatment(const double &d, const double &mu, 
                                                       const int i, const int j, const int k,
                                                       Turbulent3DWallData &wall, 
                                                       Grid3D_Hexa_Block &grid, double &dk,
                                                       double &domega){

   double omega_viscous_layer, omega_log_layer;
   
   //Check the first cells off wall only, from 2D automatic wall testing
   // it seemed that check more than one cell close to the wall
   // did not make difference numerically. (numerical solutions).
   if(((i== grid.ICl) &&(grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
                          grid.BCtypeW[j][k] == BC_NO_SLIP  ||
                          grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                          grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                          grid.BCtypeW[j][k] == BC_ADIABATIC_WALL)) 
      ||((i==grid.ICu) &&(grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
                           grid.BCtypeE[j][k] == BC_NO_SLIP  ||
                           grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                           grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                           grid.BCtypeE[j][k] == BC_ADIABATIC_WALL))
      ||((j == grid.JCl) &&(grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeS[i][k] == BC_NO_SLIP  ||
                             grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeS[i][k] == BC_ADIABATIC_WALL))
      ||((j ==grid.JCu) &&(grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
                            grid.BCtypeN[i][k] == BC_NO_SLIP  ||
                            grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeN[i][k] == BC_ADIABATIC_WALL))  
      ||((k == grid.KCl) &&(grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
                             grid.BCtypeB[i][j] == BC_NO_SLIP  ||
                             grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeB[i][j] == BC_ADIABATIC_WALL)) 
      ||((k ==grid.KCu) &&(grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
                            grid.BCtypeT[i][j] == BC_NO_SLIP  ||
                            grid.BCtypeT[i][j]  == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeT[i][j] == BC_ADIABATIC_WALL))){

       
      /* y+ <2.5 */
      if (wall.yplus <= y_sublayer){
         domega =  d*6.0*(mu/d)/(beta*sqr(wall.ywall));
         
      }
      /* y+ > 2.5 && y+<30  apply  blending formulation */
      if(wall.yplus > y_sublayer && wall.yplus <yplus_l){
         //cout<<"--------------#2"<<endl; 
         //compute omega values at both viscous sublayer and log-layer 
         omega_viscous_layer = 6.0*(mu/d)/(beta*sqr(wall.ywall)); 
         omega_log_layer = wall.utau/(sqrt(beta_star)*Karman_const*wall.ywall);
                
         //blending two omega values
         domega = d*sqrt(sqr(omega_viscous_layer)+sqr(omega_log_layer));
         //k behaves quadratically in the intermediate layer in between
         // viscous sublayer and log layer
         dk = d* sqr(wall.utau)/sqrt(beta_star)*sqr(wall.yplus/yplus_l);
                 
      }
      /* y+ >= 30  apply wall function */
      if(wall.yplus >= yplus_l) {      
         // Set k
         dk = d*sqr(wall.utau)/sqrt(beta_star);
         // Set Omega
         if(wall.utau == ZERO){
            domega = d* 6.0*(mu/d)/(beta*sqr(wall.ywall));
         }else{   
            domega = d*wall.utau/(sqrt(beta_star)*Karman_const*wall.ywall);   
            
         }
         
                  
      } /* endif */
   
   
   } //end of wall treatment switch checking 
    
   
         

   return 0;
   
   
}

int Turbulence_Model_k_omega::automatic_wall_treatment_residual_evaluaton(
   const int i, const int j, const int k,
   Turbulent3DWallData &wall,
   Grid3D_Hexa_Block &grid,
   double &dUdt_rhok,
   double &dUdt_rhoomega){


   if(((i== grid.ICl) &&(grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
                          grid.BCtypeW[j][k] == BC_NO_SLIP  ||
                          grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                          grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                          grid.BCtypeW[j][k] == BC_ADIABATIC_WALL)) 
      ||((i==grid.ICu) &&(grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
                           grid.BCtypeE[j][k] == BC_NO_SLIP  ||
                           grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                           grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                           grid.BCtypeE[j][k] == BC_ADIABATIC_WALL))
      ||((j == grid.JCl) &&(grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
                             grid.BCtypeS[i][k] == BC_NO_SLIP  ||
                             grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeS[i][k] == BC_ADIABATIC_WALL))
      ||((j ==grid.JCu) &&(grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
                            grid.BCtypeN[i][k] == BC_NO_SLIP  ||
                            grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeN[i][k] == BC_ADIABATIC_WALL))  
      ||((k == grid.KCl) &&(grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
                             grid.BCtypeB[i][j] == BC_NO_SLIP  ||
                             grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                             grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                             grid.BCtypeB[i][j] == BC_ADIABATIC_WALL)) 
      ||((k ==grid.KCu) &&(grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
                            grid.BCtypeT[i][j] == BC_NO_SLIP  ||
                            grid.BCtypeT[i][j]  == BC_WALL_VISCOUS_ISOTHERMAL ||
                            grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                            grid.BCtypeT[i][j] == BC_ADIABATIC_WALL))){
      /* y+ <2.5 */
      if (wall.yplus <= y_sublayer){
         dUdt_rhoomega =  0.0;
         
      }else{ // y>2.5
         dUdt_rhok = 0.0;
         dUdt_rhoomega = 0.0;
      }
         
   
   } //end of wall treatment switch checking 
   return 0;
   
  
}

   

/****************** EddyDissipationModel.h *****************************************
  This header file defines various classes for various standard turbulence/chemistry
  models.  
***********************************************************************************/

#ifndef _EDDY_DISSIPATION_MODELLING_INCLUDED 
#define _EDDY_DISSIPATION_MODELLING_INCLUDED 

/* Include required CFFC header files. */

#ifndef _TURBULENCE_MODELLING_INCLUDED 
#include "../TurbulenceModelling.h"
#endif //_TURBULENCE_MODELLING_INCLUDED 

#ifndef  _REACTIONS_INCLUDED
#include "../Reactions/Reactions.h"
#endif // _REACTIONS_INCLUDED 

// modeling chemistry source term with using
// eddy dissipation model for describing turbulence/chemistry interaction.
template<class SOLN_pSTATE, class SOLN_cSTATE>
SOLN_cSTATE Seddydissipationmodel(SOLN_pSTATE &W);
   
/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using eddy dissipation model.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations (c_i*rho/M_i)
              mol/m^3.

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
SOLN_cSTATE Seddydissipationmodel(SOLN_pSTATE &W) {
   
   SOLN_cSTATE U_return;
   U_return.Vacuum();

   double Temp = W.T();  //K
   double Press = W.p;    // [Pa]
   double PressDyne = W.p*TEN; // N/m^2 -> dyne/cm^2
   double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3 
   double a,b, dcdt;
   
   double model_constant = 4.0;
   double tau_c, tau_t, s ;
   
   switch(W.React.reactset_flag){
     //--------- ONE STEP CH4 ----------//
     case CH4_1STEP: 
       s = TWO*(W.specdata[1].Mol_mass()*THOUSAND)/(W.specdata[0].Mol_mass()*THOUSAND);
       if (W.spec[0].c > ZERO && W.spec[1].c > ZERO) {
          // tau_t --- s/(mol/cm^3)
          tau_t = (W.specdata[0].Mol_mass()*THOUSAND)/(W.k_omega_model.beta_star*max(TOLER, W.omega))/
                  (model_constant*rho*min(W.spec[0].c, W.spec[1].c/s) );
          tau_c +=tau_t;
       } /* endif */
     
       // compare two time scales
       //   if(kf[0]!=ZERO){
       //     tau_l =(W.rho*W.spec[0].c)/(M[0]*kf[0]*THOUSAND);
       //     tau_t = ONE/(W.beta_star*max(TOLER, W.omega)*cm1);
       //     cout<<"\n t_l = "<<tau_l<<"  t_t "<<tau_t<<endl;
       //     }
       if (tau_c > ZERO) {
          for (int index =0; index < W.ns; index++){
             switch(index) {
             case 0 : //CH4
                dcdt = -ONE/tau_c;
                break;
             case 1 : //O2
                dcdt = - TWO/tau_c;
                break;
             case 2 : //CO2
                dcdt = ONE/tau_c;
                break;
             case 3 : //H2O
                dcdt = TWO/tau_c;
                break;
             } /* endswitch */
           
             // "dcdt" in mol/(cm^3*s)
             // U_return.rhospec[].c should be in (kg/(m^3*s))
             // conversion :  g/mol * mol/(cm^3*s) = g/(10^(-6)m^3s) = 1000* kg/(m^3*s)
             U_return.rhospec[index].c =  W.specdata[index].Mol_mass()*THOUSAND*dcdt*THOUSAND;
          } /* endfor */
       } else {
          for(int index =0; index<W.ns; index++){
             switch(index) {
             case 0 : //CH4
                dcdt = ZERO;
                break;
             case 1 : //O2
                dcdt = ZERO;
                break;
             case 2 : //CO2
                dcdt = ZERO;
                break;
             case 3 : //H2O
                dcdt = ZERO;
                break;
             } /* endswitch */
           
             //dcdt in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3
             U_return.rhospec[index].c = W.specdata[index].Mol_mass()*THOUSAND*dcdt*THOUSAND;
          } /* endfor */
       } /* endif */
       break;

     default:
     break;

   } /* endswitch */
     
   return U_return;
  
}
  
#endif // _EDDY_DISSIPATION_MODELLING_INCLUDED 

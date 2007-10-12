
#ifndef _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#include "Euler3DThermallyPerfectHexaBlock.h"
#endif // _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block for a variety of multi-stage explicit *
 * time integration schemes.                            *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<Euler3D_ThermallyPerfect_pState, 
               Euler3D_ThermallyPerfect_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<Euler3D_ThermallyPerfect_pState,
                                                     Euler3D_ThermallyPerfect_cState> &IPs){

   int k_residual;
   double omega;
   
   int num_vars = NumVar();

   /* Perform update of solution variables for stage 
      i_stage of an N stage scheme. */

   /* Evaluate the time step fraction and residual storage 
      location for the stage. */
  
   switch(IPs.i_Time_Integration) {
   case TIME_STEPPING_EXPLICIT_EULER :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
        if (IPs.N_Stage == 4) { 
          if (i_stage == 4) { 
             k_residual = 0;
          } else { 
             k_residual = i_stage - 1;
          } // endif  
       } // endif 
      break;
   case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
      omega = MultiStage_Optimally_Smoothing(i_stage, 
                                             IPs.N_Stage,
                                             IPs.i_Limiter);
      k_residual = 0;
      break;
   default:
      omega = Runge_Kutta(i_stage, IPs.N_Stage);
      k_residual = 0;
      break;
   } /* endswitch */
    
   /* Update solution variables for this stage. */

   for (int k = KCl ; k <= KCu ; ++k ) {
      for (int j = JCl ; j <= JCu ; ++j ) {
         for (int i = ICl ; i <= ICu ; ++i ) {
            if (IPs.Local_Time_Stepping == 
                GLOBAL_TIME_STEPPING || 
                IPs.Local_Time_Stepping == 
                SCALAR_LOCAL_TIME_STEPPING) {
               U[i][j][k] = Uo[i][j][k] + omega* dUdt[i][j][k][k_residual];
               //N-1 species
               U[i][j][k][num_vars] = U[i][j][k].rho*(ONE - U[i][j][k].sum_species());
            } /* endif */
            
            Uo[i][j][k].negative_speccheck(Uo[i][j][k], W[i][j][k].React.reactset_flag);

            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING && 
                (U[i][j][k].rho  <= ZERO ||  
                 U[i][j][k].es() <= ZERO)) {
               cout << "\n " << CFFC_Name() 
                    << " ERROR: Negative Density, Mass Fractions, and/or Sensible Energy: \n"
                    << " cell = (" << i << ", " << j <<", "<< k << ") " 
                    << " X = " <<  Grid.Cell[i][j][k].Xc 
                    << "\n U = " <<  U[i][j][k] 
                    << "\n dUdt = " << dUdt[i][j][k][k_residual] 
                    << " omega = " << omega << "\n";
               return (i);
            } /* endif */
                          
            W[i][j][k] = U[i][j][k].W();
	 } /* endfor */    	 
      } /* endfor */    
   } /* endfor */

   /* Solution successfully updated. */
    
   return (0);   
      
}

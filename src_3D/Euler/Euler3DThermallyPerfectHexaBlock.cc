
#ifndef _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#include "Euler3DThermallyPerfectHexaBlock.h"
#endif // _EULER3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<Euler3D_ThermallyPerfect_pState, 
                Euler3D_ThermallyPerfect_cState>::
Output_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                                Euler3D_ThermallyPerfect_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {

   Euler3D_ThermallyPerfect_pState W_node;
   
   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output nodal solution data. */
   
   Out_File << setprecision(14);

   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";

      //n species mass fractions names
      for (int i =0 ; i < W[0][0][0].ns ; i++) {
         Out_File << "\"c_" << W[0][0][0].specdata[i].Speciesname()
                  << "\" \\ \n";
      } /* endfor */
      
      Out_File << "\"T\" \\ \n"
               << "\"M\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu -  Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu -  Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu -  Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   for (int k  =  Grid.KNl ; k <=  Grid.KNu ; ++k ) {
      for (int j  =  Grid.JNl ; j <=  Grid.JNu ; ++j ) {
         for (int i =  Grid.INl ; i <=  Grid.INu ; ++i ) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T() 
                     << " " << W_node.M() << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<Euler3D_ThermallyPerfect_pState, 
                Euler3D_ThermallyPerfect_cState>::
Output_Cells_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                                      Euler3D_ThermallyPerfect_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);

   /* Output cell centred solution data. */

   Out_File << setprecision(14);

   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " 
               << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";

      //n species mass fractions names
      for (int i =0; i < W[0][0][0].ns; i++) {
         Out_File << "\"c" << W[0][0][0].specdata[i].Speciesname()
                  << "\" \\ \n";
      } /* endfor */
     
      Out_File << "\"T\" \\ \n"
               << "\"M\" \\ \n";

      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  ICu -  ICl + 2*Nghost +1 << " \\ \n"
               << "J = " <<  JCu -  JCl + 2*Nghost +1 << " \\ \n"
               << "K = " <<  KCu -  KCl + 2*Nghost +1 << " \\ \n"
               << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " <<  Grid.ICu -  Grid.ICl + 2* Nghost + 1 << " \\ \n"
               << "J = " <<  Grid.JCu -  Grid.JCl + 2* Nghost + 1 << " \\ \n"
               << "K = " <<  Grid.KCu -  Grid.KCl + 2* Nghost + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";
      
   } /* endif */

   for (int k =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
      for (int j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         for (int i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            Out_File << " "  << Grid.Cell[i][j][k].Xc
                     << W[i][j][k];
            Out_File.setf(ios::scientific);
            Out_File << " " << W[i][j][k].T() 
                     << " " << W[i][j][k].M() << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
    
}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the solution values at the nodes of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template< >
void Hexa_Block<Euler3D_ThermallyPerfect_pState, 
                Euler3D_ThermallyPerfect_cState>::
Output_Nodes_Tecplot(Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                                      Euler3D_ThermallyPerfect_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {

   Euler3D_ThermallyPerfect_pState W_node;
   
   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output nodal solution data. */
   
   Out_File << setprecision(14);

   if (Output_Title) {
      Out_File << "TITLE = \"" << CFFC_Name() << ": 3D Solution, "
               << "Time Step/Iteration Level = " << Number_of_Time_Steps
               << ", Time = " << Time
               << "\"" << "\n"
               << "VARIABLES = \"x\" \\ \n"
               << "\"y\" \\ \n"
               << "\"z\" \\ \n"
               << "\"rho\" \\ \n"
               << "\"u\" \\ \n"
               << "\"v\" \\ \n"
               << "\"w\" \\ \n"
               << "\"p\" \\ \n";

      //n species mass fractions names
      for (int i =0 ; i < W[0][0][0].ns; i++) {
         Out_File << "\"c_" << W[0][0][0].specdata[i].Speciesname()
                  << "\" \\ \n";
      } /* endfor */
      
      Out_File << "\"T\" \\ \n"
               << "\"M\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu -  Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu -  Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu -  Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   for (int k  =  Grid.KNl ; k <=  Grid.KNu ; ++k ) {
      for (int j  =  Grid.JNl ; j <=  Grid.JNu ; ++j ) {
         for (int i =  Grid.INl ; i <=  Grid.INu ; ++i ) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T()
                     << " " << W_node.M() << "\n";
            Out_File.unsetf(ios::scientific);
         } /* endfor */
      } /* endfor */
   } /* endfor */

   Out_File << setprecision(6);
    
}

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

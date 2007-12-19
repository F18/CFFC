
#ifndef _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED
#include "FANS3DThermallyPerfectHexaBlock.h"
#endif // _FANS3D_THERMALLYPERFECT_HEXA_BLOCK_INCLUDED

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
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                                FANS3D_ThermallyPerfect_KOmega_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {

   FANS3D_ThermallyPerfect_KOmega_pState W_node;
   
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
               << "\"p\" \\ \n"
               << "\"k\" \\ \n"
               << "\"omega\" \\ \n";

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
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values of the       *
 * specified hexadedral solution block to the           *
 * specified output stream suitable for plotting with   *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Cells_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs,
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
               << "\"p\" \\ \n"
               << "\"k\" \\ \n"
               << "\"omega\" \\ \n";

      //n species mass fractions names
      for (int i = 0; i < W[0][0][0].ns; i++) {
         Out_File << "\"c" << W[0][0][0].specdata[i].Speciesname()
                  << "\" \\ \n";
      } /* endfor */
     
      Out_File <<"\"T\" \\ \n"
               <<"\"M\" \\ \n"
               << "\"ywall\" \\ \n"
               << "\"yplus\" \\ \n"
               << "\"tauw\" \\ \n"
               << "\" utau\" \\ \n";
      
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
            Out_File << " " << Grid.Cell[i][j][k].Xc
                     << W[i][j][k];
            Out_File.setf(ios::scientific);
            Out_File << " " << W[i][j][k].T() 
                     << " " << W[i][j][k].M();
            Out_File << " " << WallData[i][j][k].ywall 
                     << " " << WallData[i][j][k].yplus
                     << " " << WallData[i][j][k].tauw 
                     << " " << WallData[i][j][k].utau << "\n ";
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
template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
Output_Nodes_Tecplot(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                                      FANS3D_ThermallyPerfect_KOmega_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {

   FANS3D_ThermallyPerfect_KOmega_pState W_node;
   
   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);
   
   /* Output nodal solution data. */
   
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
               << "\"p\" \\ \n"
               << "\"k\" \\ \n"
               << "\"omega\" \\ \n";

      //n species mass fractions names
      for (int i = 0; i < W[0][0][0].ns; i++) {
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
 * Routine: ICs                                         *
 *                                                      *
 * Apply initial conditions for the specified hexa      *
 * solution block.                                      *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
ICs(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs) {
   
   double dpdx, dpdy, dpdz, delta_pres, delta_pres_x, delta_pres_y, delta_pres_z;
   double zd, zz, di, Um, U_axi;
   
   double Rprime, yprime, xn, yn, fc, tempvalue;
   Vector2D Xt;
   ReactiveScalarField_Fuel_CH4H2 RSF;
   NonreactiveScalarField NRSF;
   double tc[10], max_mean_velocity, 
      BluffBody_Coflow_Air_Velocity = 125, 
      BluffBody_Coflow_Fuel_Velocity = 25;
   
   int BluffBody_Data_Usage = 0;

   FANS3D_ThermallyPerfect_KOmega_pState Wl, Wr;

   switch(IPs.i_ICs) {
      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X :
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.Grid_IP.Box_Width;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.x = HALF*(-dpdx)/W[i][j][k].mu()*
                                   (Grid.Cell[i][j][k].Xc.y + 0.5*IPs.Grid_IP.Box_Height)*
                                   (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.Grid_IP.Box_Height) +
                                   IPs.Moving_Wall_Velocity.x*
                                   (Grid.Cell[i][j][k].Xc.y/IPs.Grid_IP.Box_Height+HALF);
                  W[i][j][k].p = IPs.Wo.p - (i-ICl-Nghost)*delta_pres/(ICu-ICl);	 
                  if (i == ICl-Nghost || i == ICl-Nghost+1 ) {
                    W[i][j][k].p = IPs.Wo.p;
                  } /* endif */
                  if (i == ICu+Nghost || i == ICu+Nghost-1 ) {
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  } /* endif */
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y :
         dpdy = IPs.Pressure_Gradient.y;  
         delta_pres = dpdy*IPs.Grid_IP.Box_Height;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.y = HALF*(-dpdy)/W[i][j][k].mu()*
                                   (Grid.Cell[i][j][k].Xc.x + 0.5*IPs.Grid_IP.Box_Width)*
                                   (Grid.Cell[i][j][k].Xc.x - 0.5*IPs.Grid_IP.Box_Width) +
                                   IPs.Moving_Wall_Velocity.y*
                                   (Grid.Cell[i][j][k].Xc.x/IPs.Grid_IP.Box_Width+HALF);
                  W[i][j][k].p = IPs.Wo.p - (j-JCl-Nghost)*delta_pres/(JCu-JCl);
                  if (j == JCl-Nghost || j == JCl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  } /* endif */
                  if (j == JCu+Nghost || j == JCu+Nghost-1 ){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  } /* endif */
                  U[i][j][k] = W[i][j][k].U();
 	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 

      case IC_VISCOUS_COUETTE :
      case IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z :
         dpdz = IPs.Pressure_Gradient.z;  
         delta_pres = dpdz*IPs.Grid_IP.Box_Length;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.z = HALF*(-dpdz)/W[i][j][k].mu()*
                                   (Grid.Cell[i][j][k].Xc.y + 0.5*IPs.Grid_IP.Box_Height)*
                                   (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.Grid_IP.Box_Height) +
                                   IPs.Moving_Wall_Velocity.z*
                                   (Grid.Cell[i][j][k].Xc.y/IPs.Grid_IP.Box_Height+HALF);
                  W[i][j][k].p = IPs.Wo.p - (k-KCl-Nghost)*delta_pres/(KCu-KCl);
                  if (k == KCl-Nghost || k == KCl-Nghost+1 ) {
                     W[i][j][k].p = IPs.Wo.p;
                  } /* endif */
                  if (k == KCu+Nghost || k == KCu+Nghost-1 ) {
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  } /* endif */
                  U[i][j][k] = W[i][j][k].U( );
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break; 
      
      case IC_PRESSURE_GRADIENT_X :
         dpdx = IPs.Pressure_Gradient.x;  
         delta_pres = dpdx*IPs.Grid_IP.Box_Width;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.x = -HALF*(dpdx/W[i][j][k].mu())*
                                   (Grid.Cell[i][j][k].Xc.y + 0.5*IPs.Grid_IP.Box_Height)*
                                   (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.Grid_IP.Box_Height);
                  W[i][j][k].p = IPs.Wo.p - (i-ICl-Nghost)*delta_pres/(ICu-ICl);	 
                  if( i == ICl-Nghost || i == ICl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  } /* endif */
                  if( i == ICu+Nghost || i == ICu+Nghost-1){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  } /* endif */
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */ 
	    } /* endfor */       
	 } /* endfor */
         break; 

      case IC_PRESSURE_GRADIENT_Y :
         dpdx = IPs.Pressure_Gradient.y;  
         delta_pres = dpdx*IPs.Grid_IP.Box_Height;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.y = -HALF*(dpdx/W[i][j][k].mu())*
                                   (Grid.Cell[i][j][k].Xc.x + 0.5*IPs.Grid_IP.Box_Width)*
                                   (Grid.Cell[i][j][k].Xc.x - 0.5*IPs.Grid_IP.Box_Width);
                  W[i][j][k].p = IPs.Wo.p - (j-JCl-Nghost)*delta_pres/(JCu-JCl);	 
                  if( j == JCl-Nghost || j == JCl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  if( j == JCu+Nghost || j == JCu+Nghost-1){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */ 
	    } /* endfor */       
	 } /* endfor */
         break;

      case IC_PRESSURE_GRADIENT_Z :
         dpdx = IPs.Pressure_Gradient.z;  
         delta_pres = dpdx*IPs.Grid_IP.Box_Length;
         for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
            for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.z = -HALF*(dpdx/W[i][j][k].mu())*
                                   (Grid.Cell[i][j][k].Xc.y + 0.5*IPs.Grid_IP.Box_Height)*
                                   (Grid.Cell[i][j][k].Xc.y - 0.5*IPs.Grid_IP.Box_Height);
                  W[i][j][k].p = IPs.Wo.p - (k-KCl-Nghost)*delta_pres/(KCu-KCl);	 
                  if( k == KCl-Nghost || k == KCl-Nghost+1 ){
                     W[i][j][k].p = IPs.Wo.p;
                  }
                  if( k == KCu+Nghost || k == KCu+Nghost-1){
                     W[i][j][k].p = IPs.Wo.p - delta_pres; 
                  }
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */ 
	    } /* endfor */       
	 } /* endfor */
         break;

      case IC_CHANNEL_FLOW :
         // John Laufer, NACA case: 
         // Investigation of Turbulent Flow in a Two-Dimensional Channel
         dpdx = IPs.Pressure_Gradient.x;  
         dpdy = IPs.Pressure_Gradient.y;  
         dpdz = IPs.Pressure_Gradient.z;  
         delta_pres_x = dpdx*IPs.Grid_IP.Box_Width;
         delta_pres_y = dpdy*IPs.Grid_IP.Box_Height;
         delta_pres_z = dpdz*IPs.Grid_IP.Box_Length;
         Um = IPs.Reynolds_Number*IPs.Wo.mu()/(IPs.Wo.rho*IPs.Grid_IP.Box_Height);
         if (IPs.Pressure_Gradient.x != ZERO) {
            for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	       for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
                  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                     W[i][j][k] = IPs.Wo;
                     WallData[i][j][k].tauw = fabs(-IPs.Grid_IP.Box_Height*dpdx);
                     WallData[i][j][k].utau = sqrt(WallData[i][j][k].tauw/W[i][j][k].rho);
                     W[i][j][k].k = WallData[i][j][k].utau*WallData[i][j][k].utau/
                                    sqrt(W[0][0][0].k_omega_model.beta_star);
                     W[i][j][k].p = IPs.Wo.p - 
                                    (Grid.Cell[i][j][k].Xc.x)*delta_pres_x/IPs.Grid_IP.Box_Width;	 
                     // setting the u velocity in a turbulent channel flow
                     // by using the power law of u velocity in a turbulent pipe flow
                     zd = Grid.Cell[i][j][k].Xc.y;
                     if (zd > 0.0 && zd < 0.5*IPs.Grid_IP.Box_Height) {
                         zz = (1.0- zd/(0.5*IPs.Grid_IP.Box_Height));
                        W[i][j][k].v.x = Um*pow(zz, 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = ( 0.5886 - 31.2*(zd)+2039.6*zd*zd -19208*zd*zd*zd );
                        if(zd <=0.01) {W[i][j][k].k = 0.465;}
                     } else if(zd < ZERO && fabs(zd)< 0.5*IPs.Grid_IP.Box_Height) {
                        zz = (fabs(zd)/(0.5*IPs.Grid_IP.Box_Height) - 1.0);
                        W[i][j][k].v.x = Um*pow(fabs(zz), 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = (0.5886 - 31.2*fabs(zd)+2039.6*zd*zd -
                                       19208*fabs(zd)*zd*zd );
                        if(fabs(zd) <=0.01) {W[i][j][k].k = 0.465;}
                     } else {
                        W[i][j][k].v.x = Um;
                     } /* endif */
                     if (WallData[i][j][k].ywall !=0.0) {
                        W[i][j][k].omega = WallData[i][j][k].utau/
                                           (sqrt(W[0][0][0].k_omega_model.beta_star)*
                                           W[0][0][0].k_omega_model.Karman_const*
                                           WallData[i][j][k].ywall); 
		     } /* endif */
                     //conservative solution state
                     U[i][j][k] = W[i][j][k].U();
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */

         } else if (IPs.Pressure_Gradient.y !=ZERO) {
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
	       for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
                  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                     W[i][j][k] = IPs.Wo;
                     WallData[i][j][k].tauw = fabs(-IPs.Grid_IP.Box_Width*dpdy);
                     WallData[i][j][k].utau = sqrt(WallData[i][j][k].tauw/W[i][j][k].rho);
                     W[i][j][k].k = WallData[i][j][k].utau*WallData[i][j][k].utau/
                                    sqrt(W[0][0][0].k_omega_model.beta_star);
                     W[i][j][k].p = IPs.Wo.p - 
                                    (Grid.Cell[i][j][k].Xc.y)*delta_pres_y/IPs.Grid_IP.Box_Height;	 
                     // setting the u velocity in a turbulent channel flow
                     // by using the power law of u velocity in a turbulent pipe flow
                     zd = Grid.Cell[i][j][k].Xc.x;
                     if (zd > 0.0 && zd < 0.5*IPs.Grid_IP.Box_Width) {
                        zz = (1.0- zd/(0.5*IPs.Grid_IP.Box_Width));
                        W[i][j][k].v.y = Um*pow(zz, 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = ( 0.5886 - 31.2*(zd)+2039.6*zd*zd -19208*zd*zd*zd );
                        if(zd <=0.01) {W[i][j][k].k = 0.465;}
                     } else if (zd<ZERO && fabs(zd)< 0.5*IPs.Grid_IP.Box_Width) {
                        zz = (fabs(zd)/(0.5*IPs.Grid_IP.Box_Width) - 1.0);
                        W[i][j][k].v.y = Um*pow(fabs(zz), 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = (0.5886 - 31.2*fabs(zd)+2039.6*zd*zd -
                                        19208*fabs(zd)*zd*zd );
                        if (fabs(zd) <=0.01) { W[i][j][k].k = 0.465;}
                     } else {
                        W[i][j][k].v.y = Um;
		     } /* endif */
                     if(WallData[i][j][k].ywall !=0.0){
                        W[i][j][k].omega = WallData[i][j][k].utau/
                                           (sqrt(W[0][0][0].k_omega_model.beta_star)*
                                           W[0][0][0].k_omega_model.Karman_const*
                                           WallData[i][j][k].ywall); 
		     } /* endif */
                     //conservative solution state
                     U[i][j][k] = W[i][j][k].U();
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */

         } else if (IPs.Pressure_Gradient.z !=ZERO) {
            for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
	       for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
                  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                     W[i][j][k] = IPs.Wo;
                     WallData[i][j][k].tauw = fabs(-IPs.Grid_IP.Box_Height*dpdz);
                     WallData[i][j][k].utau = sqrt(WallData[i][j][k].tauw/W[i][j][k].rho);
                     W[i][j][k].k = WallData[i][j][k].utau*WallData[i][j][k].utau/
                                    sqrt(W[0][0][0].k_omega_model.beta_star);
                     W[i][j][k].p = IPs.Wo.p - 
                                    (Grid.Cell[i][j][k].Xc.z)*delta_pres_z/IPs.Grid_IP.Box_Length;	 
                     // setting the u velocity in a turbulent channel flow
                     // by using the power law of u velocity in a turbulent pipe flow
                     zd = Grid.Cell[i][j][k].Xc.y;
                     if (zd > 0.0 && zd < 0.5*IPs.Grid_IP.Box_Height) {
                         zz = (1.0- zd/(0.5*IPs.Grid_IP.Box_Height));
                        W[i][j][k].v.z = Um*pow(zz, 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = (0.5886 - 31.2*(zd)+2039.6*zd*zd -19208*zd*zd*zd);
                        if (zd <=0.01) {W[i][j][k].k = 0.465;}
                     } else if (zd<ZERO && fabs(zd)< 0.5*IPs.Grid_IP.Box_Height) {
                        zz = (fabs(zd)/(0.5*IPs.Grid_IP.Box_Height) - 1.0);
                        W[i][j][k].v.z = Um*pow(fabs(zz), 0.133);
                        // try to feed a reasonalbe k profile as initial condition
                        // to see how far it goes, since this case has really small 
                        // pressure gradient
                        W[i][j][k].k = ( 0.5886 - 31.2*fabs(zd)+2039.6*zd*zd -
                                       19208*fabs(zd)*zd*zd );
                        if (fabs(zd) <= 0.01) {W[i][j][k].k = 0.465;}
                     } else {
                        W[i][j][k].v.z = Um;
                     } /* endif */
                     if (WallData[i][j][k].ywall !=0.0) {
                        W[i][j][k].omega = WallData[i][j][k].utau/
                                           (sqrt(W[0][0][0].k_omega_model.beta_star)*
                                           W[0][0][0].k_omega_model.Karman_const*
                                           WallData[i][j][k].ywall); 
		     } /* endif */
                     //conservative solution state
                     U[i][j][k] = W[i][j][k].U();
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */

	 } /* endif */         
         break; 

      case IC_SHOCK_BOX :
	 Wl = FANS3D_ThermallyPerfect_KOmega_pState(IPs.Wo);
         Wl.rho = DENSITY_STDATM;
         Wl.v = Vector3D_ZERO;
         Wl.p = PRESSURE_STDATM;
         Wr = FANS3D_ThermallyPerfect_KOmega_pState(IPs.Wo);
         Wr.rho = EIGHT*DENSITY_STDATM;
         Wr.v = Vector3D_ZERO;
         Wr.p = TEN*PRESSURE_STDATM; 
         for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
	    for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
               for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
		  if (Grid.Cell[i][j][k].Xc.x <= ZERO &&
                      Grid.Cell[i][j][k].Xc.y <= ZERO &&
                      Grid.Cell[i][j][k].Xc.z <= ZERO) {
                      W[i][j][k] = Wl; 
                  } else {
                      W[i][j][k] = Wr;
                  } /* end if */
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break;

      case IC_SOD_XDIR :
         Wl = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM, 
                                                    Vector3D_ZERO,
                                                    PRESSURE_STDATM,
                                                    ZERO,
                                                    ZERO);
         Wr = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM*8.0, 
                                                    Vector3D_ZERO,
                                                    PRESSURE_STDATM*10.0,
                                                    ZERO,
                                                    ZERO);
	 for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
            for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
               for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
                  if (Grid.Cell[i][j][k].Xc.x <= ZERO) {
                     W[i][j][k] = Wl; 
                  } else {
                     W[i][j][k] = Wr;
		  } /* end if */
                  U[i][j][k] = W[i][j][k].U();
               } /* endfor */
            } /* endfor */
         } /* endfor */
	 break;

      case IC_SOD_YDIR :
         Wl = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM, 
                                                    Vector3D_ZERO,
	                                            PRESSURE_STDATM,
                                                    ZERO,
                                                    ZERO);
         Wr = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM*8.0, 
                                                    Vector3D_ZERO,
                                                    PRESSURE_STDATM*10.0,
                                                    ZERO,
                                                    ZERO);
         for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
            for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
               for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
                  if (Grid.Cell[i][j][k].Xc.y <= ZERO) {
                     W[i][j][k] = Wl; 
		  } else {
                     W[i][j][k] = Wr;
                  } /* end if */
                  U[i][j][k] = W[i][j][k].U();
               } /* endfor */
            } /* endfor */
         } /* endfor */
         break;

      case IC_SOD_ZDIR :
         Wl = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM, 
                                                    Vector3D_ZERO,
                                                    PRESSURE_STDATM,
                                                    ZERO,
                                                    ZERO);
         Wr = FANS3D_ThermallyPerfect_KOmega_pState(DENSITY_STDATM*8.0, 
                                                    Vector3D_ZERO,
                                                    PRESSURE_STDATM*10.0,
                                                    ZERO,
                                                    ZERO);
         for (int k = KCl-Nghost; k <= KCu+Nghost; ++k) {
            for (int j = JCl-Nghost; j <= JCu+Nghost; ++j) {
               for (int i = ICl-Nghost; i <= ICu+Nghost; ++i) {
                  if (Grid.Cell[i][j][k].Xc.z <= ZERO) {
                     W[i][j][k] = Wl; 
                  } else {
                     W[i][j][k] = Wr;
		  } /* end if */
                  U[i][j][k] = W[i][j][k].U();
               } /* endfor */
            } /* endfor */
         } /* endfor */
         break;

      case IC_TURBULENT_COFLOW:
         for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
            for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                  // Apply uniform solution state
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.z= ZERO;
                  W[i][j][k].v.x = ZERO;
                  W[i][j][k].v.y = ZERO;
       
                  W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300.0);
                  tempvalue =  BluffBody_Coflow_Air_Velocity*(HALF*(IPs.Grid_IP.Radius_Coflow_Inlet_Pipe-
                               IPs.Grid_IP.Radius_Bluff_Body))/(W[i][j][k].mu()/W[i][j][k].rho);
                  WallData[i][j][k].tauw = 0.0228*W[i][j][k].rho* BluffBody_Coflow_Air_Velocity*
                                           BluffBody_Coflow_Air_Velocity/pow(fabs(tempvalue), 0.25);
                  WallData[i][j][k].utau = sqrt(WallData[i][j][k].tauw/W[i][j][k].rho);
                  W[i][j][k].k = 0.5*WallData[i][j][k].utau*WallData[i][j][k].utau/
                                 sqrt(W[0][0][0].k_omega_model.beta_star);
                  WallData[i][j][k].yplus = WallData[i][j][k].ywall*WallData[i][j][k].utau/
                                            (W[i][j][k].mu()/W[i][j][k].rho);
                  if (WallData[i][j][k].ywall !=ZERO) {
                     if (WallData[i][j][k].yplus<=W[0][0][0].k_omega_model.y_sublayer) {
                        W[i][j][k].omega = 6.0*(W[i][j][k].mu()/W[i][j][k].rho)/
                                           (W[0][0][0].k_omega_model.beta*WallData[i][j][k].ywall*WallData[i][j][k].ywall);
		     } else {
                        W[i][j][k].omega = sqrt(W[i][j][k].k)/(pow(W[0][0][0].k_omega_model.beta_star, 0.25)*
                                           W[0][0][0].k_omega_model.Karman_const*WallData[i][j][k].ywall);
		     }
		  }
                  // Specifying the velocity profiles in the annular pipe (coflowing air and jet)
                  if ((fabs(Grid.Cell[i][j][k].Xc.y)>IPs.Grid_IP.Radius_Bluff_Body) && 
                      (fabs(Grid.Cell[i][j][k].Xc.y)<IPs.Grid_IP.Radius_Coflow_Inlet_Pipe)) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = (IPs.Grid_IP.Radius_Coflow_Inlet_Pipe - IPs.Grid_IP.Radius_Bluff_Body)/2.0;
                        if (fabs(Grid.Cell[i][j][k].Xc.y)<=(IPs.Grid_IP.Radius_Bluff_Body+ Rprime)) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) -  IPs.Grid_IP.Radius_Bluff_Body;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) -  IPs.Grid_IP.Radius_Bluff_Body- Rprime;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity;
                           }
			}
		     }
		  }
                  if ((fabs(Grid.Cell[i][j][k].Xc.x)>IPs.Grid_IP.Radius_Bluff_Body) && 
                      (fabs(Grid.Cell[i][j][k].Xc.x)<IPs.Grid_IP.Radius_Coflow_Inlet_Pipe)) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = (IPs.Grid_IP.Radius_Coflow_Inlet_Pipe - IPs.Grid_IP.Radius_Bluff_Body)/2.0;
                        if (fabs(Grid.Cell[i][j][k].Xc.x)<=(IPs.Grid_IP.Radius_Bluff_Body+ Rprime)) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) -  IPs.Grid_IP.Radius_Bluff_Body;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity;
			   }
			} else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) -  IPs.Grid_IP.Radius_Bluff_Body- Rprime;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Air_Velocity;
			   }
			}
		     }
		  }
                  if (fabs(Grid.Cell[i][j][k].Xc.y) <= IPs.Grid_IP.Radius_Fuel_Line) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = IPs.Grid_IP.Radius_Fuel_Line;
                        if (Grid.Cell[i][j][k].Xc.y<0.0) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y);
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) ;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                           }
                        }
                     }
		  }
                  if (fabs(Grid.Cell[i][j][k].Xc.x) <= IPs.Grid_IP.Radius_Fuel_Line) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) { 
                        Rprime = IPs.Grid_IP.Radius_Fuel_Line;
                        if (Grid.Cell[i][j][k].Xc.x<0.0) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x);
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) ;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                           }
                        }
		     }
		  }
                  if (IPs.Species_IP.num_species == 3) { 
                     if (fabs(Grid.Cell[i][j][k].Xc.x)<0.5*IPs.Grid_IP.Radius_Fuel_Line && 
                         Grid.Cell[i][j][k].Xc.z>0.0 && 
                         Grid.Cell[i][j][k].Xc.z<0.1*IPs.Grid_IP.Length_Combustor_Tube) {
                        W[i][j][k].spec[0] = 1.0; //CH4
                        W[i][j][k].spec[1] = 0.0;//O2
                        W[i][j][k].spec[2] = 0.0; //N2
                     }
                     if (fabs(Grid.Cell[i][j][k].Xc.y)<0.5*IPs.Grid_IP.Radius_Fuel_Line && 
                         Grid.Cell[i][j][k].Xc.z>0.0 && 
                         Grid.Cell[i][j][k].Xc.z<0.1*IPs.Grid_IP.Length_Combustor_Tube) {
                        W[i][j][k].spec[0] = 1.0; //CH4
                        W[i][j][k].spec[1] = 0.0;//O2
                        W[i][j][k].spec[2] = 0.0; //N2
                     }
		  }
                  W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300.0);
                  U[i][j][k] = W[i][j][k].U();
	       } /* endfor */
	    } /* endfor */
         } /* endfor */
         break;

      case IC_TURBULENT_DIFFUSION_FLAME :  
         for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
            for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
                  // Apply uniform solution state
                  W[i][j][k] = IPs.Wo;
                  W[i][j][k].v.z= ZERO;
                  W[i][j][k].v.x = ZERO;
                  W[i][j][k].v.y = ZERO;
                  W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300.0);
                  tempvalue =  BluffBody_Coflow_Air_Velocity*(HALF*(IPs.Grid_IP.Radius_Coflow_Inlet_Pipe-
                               IPs.Grid_IP.Radius_Bluff_Body))/(W[i][j][k].mu()/W[i][j][k].rho);
                  WallData[i][j][k].tauw = 0.0228*W[i][j][k].rho* BluffBody_Coflow_Air_Velocity*
                                           BluffBody_Coflow_Air_Velocity/pow(fabs(tempvalue), 0.25);
                  WallData[i][j][k].utau = sqrt(WallData[i][j][k].tauw/W[i][j][k].rho);
                  W[i][j][k].k = 0.5*WallData[i][j][k].utau*WallData[i][j][k].utau/sqrt(W[0][0][0].k_omega_model.beta_star);
                  WallData[i][j][k].yplus = WallData[i][j][k].ywall*WallData[i][j][k].utau/
                                            (W[i][j][k].mu()/W[i][j][k].rho);
                  if (WallData[i][j][k].ywall !=ZERO) {
                     if (WallData[i][j][k].yplus<=W[0][0][0].k_omega_model.y_sublayer) {
                        W[i][j][k].omega = 6.0*(W[i][j][k].mu()/W[i][j][k].rho)/
                                           (W[0][0][0].k_omega_model.beta*WallData[i][j][k].ywall*WallData[i][j][k].ywall);
                     } else {
                        W[i][j][k].omega = sqrt(W[i][j][k].k)/(pow(W[0][0][0].k_omega_model.beta_star, 0.25)*
                                           W[0][0][0].k_omega_model.Karman_const*WallData[i][j][k].ywall);
                     }
                  }
                  // Specifying the velocity profiles in the annular pipe (coflowing air and jet)
                  if ((fabs(Grid.Cell[i][j][k].Xc.y)>IPs.Grid_IP.Radius_Bluff_Body) && 
                      (fabs(Grid.Cell[i][j][k].Xc.y)<IPs.Grid_IP.Radius_Coflow_Inlet_Pipe)) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = (IPs.Grid_IP.Radius_Coflow_Inlet_Pipe - IPs.Grid_IP.Radius_Bluff_Body)/2.0;
                        if (fabs(Grid.Cell[i][j][k].Xc.y)<=(IPs.Grid_IP.Radius_Bluff_Body+ Rprime)){
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) -  IPs.Grid_IP.Radius_Bluff_Body;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) -  IPs.Grid_IP.Radius_Bluff_Body- Rprime;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity;
                           }
                        }
                     }
                  }
                  if ((fabs(Grid.Cell[i][j][k].Xc.x)>IPs.Grid_IP.Radius_Bluff_Body) && 
                      (fabs(Grid.Cell[i][j][k].Xc.x)<IPs.Grid_IP.Radius_Coflow_Inlet_Pipe)) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = (IPs.Grid_IP.Radius_Coflow_Inlet_Pipe - IPs.Grid_IP.Radius_Bluff_Body)/2.0;
                        if (fabs(Grid.Cell[i][j][k].Xc.x)<=(IPs.Grid_IP.Radius_Bluff_Body+ Rprime)) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) -  IPs.Grid_IP.Radius_Bluff_Body;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) -  IPs.Grid_IP.Radius_Bluff_Body- Rprime;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Air_Velocity;
                           }
                        }
                     }
                  }
                  if (fabs(Grid.Cell[i][j][k].Xc.y) <= IPs.Grid_IP.Radius_Fuel_Line) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = IPs.Grid_IP.Radius_Fuel_Line;
                        if (Grid.Cell[i][j][k].Xc.y<0.0) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y);
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Fuel_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Fuel_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.y) ;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z = BluffBody_Coflow_Fuel_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                           } else {
                              W[i][j][k].v.z = BluffBody_Coflow_Fuel_Velocity;
                           }
                        }
                     }
                  }
                  if (fabs(Grid.Cell[i][j][k].Xc.x) <= IPs.Grid_IP.Radius_Fuel_Line) {
                     if (Grid.Cell[i][j][k].Xc.z <=  0.5*IPs.Grid_IP.Length_Combustor_Tube) {
                        Rprime = IPs.Grid_IP.Radius_Fuel_Line;
                        if (Grid.Cell[i][j][k].Xc.x<0.0) {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x);
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(yprime/Rprime), 0.143);
                           } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                           }
                        } else {
                           yprime = fabs(Grid.Cell[i][j][k].Xc.x) ;
                           if (yprime !=0.0) {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity*pow(fabs(1.0- fabs(yprime/Rprime)), 0.143);
                          } else {
                              W[i][j][k].v.z =  BluffBody_Coflow_Fuel_Velocity;
                          }
			}
                     }
                  }
                  if (IPs.Wo.React.reactset_flag == CH4_1STEP) {
                     if (fabs(Grid.Cell[i][j][k].Xc.x)<0.5*IPs.Grid_IP.Radius_Fuel_Line && 
                         Grid.Cell[i][j][k].Xc.z>0.0 && 
                         Grid.Cell[i][j][k].Xc.z<0.1*IPs.Grid_IP.Length_Combustor_Tube) {
                        double profile = 0.996;
                        W[i][j][k].spec[0] = profile;
                        W[i][j][k].spec[1] = (ONE-profile)*0.235;
                        W[i][j][k].spec[2] = ZERO; //CO2
                        W[i][j][k].spec[3] = ZERO ; //H2O
                        W[i][j][k].spec[ IPs.Wo.ns-1] =(ONE-profile)*0.765;//N2
                        W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300);
                     }
                     if (fabs(Grid.Cell[i][j][k].Xc.y)<0.5*IPs.Grid_IP.Radius_Fuel_Line && 
                         Grid.Cell[i][j][k].Xc.z>0.0 && 
                         Grid.Cell[i][j][k].Xc.z<0.1*IPs.Grid_IP.Length_Combustor_Tube) {
                        double profile = 0.996;// ONE*pow((1.0-abs(Grid.Cell[i][j][k].Xc.y)/IPs.Grid_IP.Radius_Fuel_Line), 0.143);
                        W[i][j][k].spec[0] = profile;
                        W[i][j][k].spec[1] = (ONE-profile)*0.235;
                        W[i][j][k].spec[2] = ZERO; //CO2
                        W[i][j][k].spec[3] = ZERO ; //H2O
                        W[i][j][k].spec[ IPs.Wo.ns-1] =(ONE-profile)*0.765;//N2
                        W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300);
                     }
                  }// mass fraction distribution for CH4 bluff body burner at the fuel inlet
                  W[i][j][k].rho = W[i][j][k].p/( W[i][j][k].Rtot()*300.0);
                  U[i][j][k] = W[i][j][k].U();
               } /* endfor */
	    } /* endfor */
         } /* endfor */
         break;

      case IC_UNIFORM :
      default:
         // Set the solution state everywhere to the initial state Wo[0].
	 for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	    for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
               for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
                  W[i][j][k] = IPs.Wo;
                  U[i][j][k] = W[i][j][k].U( );
               } /* endfor */
	    } /* endfor */
	 } /* endfor */
         break;
      
   } /* endswitch */

   /* Set default values for the boundary conditions
      reference states. */

   for (int k = KCl-Nghost ; k<= KCu+Nghost; ++k ) {
      for (int j = JCl-Nghost ; j<= JCu+Nghost; ++j ){
         if ((k >= KCl && k <= KCu) && (j >= JCl && j <= JCu)) {
            WoW[j][k] = W[ICl][j][k];
            WoE[j][k] = W[ICu][j][k];
         } else if (j < JCl && k < KCl ) {
            WoW[j][k] = W[ICl][JCl][KCl];
            WoE[j][k] = W[ICu][JCl][KCl];
         } else if (j > JCu && k> KCu) {
            WoW[j][k] = W[ICl][JCu][KCu];
            WoE[j][k] = W[ICu][JCu][KCu];
         } else if(j < JCl &&(k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCl][k];
            WoE[j][k] = W[ICu][JCl][k];
         } else if(j > JCu && (k >= KCl && k <= KCu)){
            WoW[j][k] = W[ICl][JCu][k];
            WoE[j][k] = W[ICu][JCu][k];
         } else if(k < KCl &&(j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCl];
            WoE[j][k] = W[ICu][j][KCl];
         } else if(k > KCu && (j >= JCl && j <= JCu)){
            WoW[j][k] = W[ICl][j][KCu];
            WoE[j][k] = W[ICu][j][KCu];
         } else if(k > KCu && j < JCl ){
            WoW[j][k] = W[ICl][JCl][KCu];
            WoE[j][k] = W[ICu][JCl][KCu];
         } else if(k < KCl && j > JCu){
            WoW[j][k] = W[ICl][JCu][KCl];
            WoE[j][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */ 
   } /* endfor */
    
   for (int k = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
         if ((k >= KCl && k <= KCu) && (i >= ICl && i <= ICu)) {
            WoS[i][k] = W[i][JCl][k];
            WoN[i][k] = W[i][JCu][k];
         } else if (i < ICl && k< KCl) {
            WoS[i][k] = W[ICl][JCl][KCl];
            WoN[i][k] = W[ICl][JCu][KCl];
         } else if (i > ICu && k > KCu) {
            WoS[i][k] = W[ICu][JCl][KCu];
            WoN[i][k] = W[ICu][JCu][KCu];
         } else if (i<ICl && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICl][JCl][k];
            WoN[i][k] = W[ICl][JCu][k];
         } else if (i>ICu && (k >= KCl && k <= KCu)){
            WoS[i][k] = W[ICu][JCl][k];
            WoN[i][k] = W[ICu][JCu][k];
         } else if ((i >= ICl && i <= ICu) && k< KCl) {
            WoS[i][k] = W[i][JCl][KCl];
            WoN[i][k] = W[i][JCu][KCl];
         } else if ((i >= ICl && i <= ICu) && k > KCu) {
            WoS[i][k] = W[i][JCl][KCu];
            WoN[i][k] = W[i][JCu][KCu];
         } else if (i < ICl  && k > KCu) {
            WoS[i][k] = W[ICl][JCl][KCu];
            WoN[i][k] = W[ICl][JCu][KCu];
         } else if (i >ICu  && k < KCl) {
            WoS[i][k] = W[ICu][JCl][KCl];
            WoN[i][k] = W[ICu][JCu][KCl];
         }
      } /* endfor */
   } /* endfor */

   for (int j = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          if ((j >= JCl && j <= JCu) && (i >= ICl && i <= ICu)) {
             WoT[i][j] = W[i][j][KCu];
             WoB[i][j] = W[i][j][KCl];
          } else if (i < ICl &&  j< JCl) {
             WoT[i][j] = W[ICl][JCl][KCu];
             WoB[i][j] = W[ICl][JCl][KCl];
          } else if(i > ICu &&  j > JCu) {
             WoT[i][j] = W[ICu][JCu][KCu];
             WoB[i][j] = W[ICu][JCu][KCl];
          }else if (i < ICl && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICl][j][KCu];
             WoB[i][j] = W[ICl][j][KCl];
          }else if (i > ICu && (j >= JCl && j <= JCu)) {
             WoT[i][j] = W[ICu][j][KCu];
             WoB[i][j] = W[ICu][j][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j< JCl) {
             WoT[i][j] = W[i][JCl][KCu];
             WoB[i][j] = W[i][JCl][KCl];
          } else if ((i >= ICl && i <= ICu) &&  j> JCu) {
             WoT[i][j] = W[i][JCu][KCu];
             WoB[i][j] = W[i][JCu][KCl];
          } else if (i > ICu && j < JCl) {
             WoT[i][j] = W[ICu][JCl][KCu];
             WoB[i][j] = W[ICu][JCl][KCl];
          } else if (i < ICl && j > JCu) {
             WoT[i][j] = W[ICl][JCu][KCu];
             WoB[i][j] = W[ICl][JCu][KCl];
          }
       } /* endfor */
   } /* endfor */
      
   // compute y+ etc.;
   Wall_Shear();

   return (0);
    
}

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified hexa solution block.                       *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                FANS3D_ThermallyPerfect_KOmega_cState>::
BCs(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs) {
   
   double dpdx, dpdy, dpdz;
   Vector3D dX;

   for (int k = KCl - Nghost; k <= KCu + Nghost; ++k) {
      for (int j = JCl - Nghost; j <= JCu + Nghost; ++j ) {
         // Prescribe West boundary conditions.
         switch(Grid.BCtypeW[j][k]) {
          case BC_NONE :
            break;
            
          case BC_REFLECTION :
            W[ICl-1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[ICl][j][k],
                                                                            Grid.nfaceW(ICl,j,k));
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[ICl+1][j][k],
                                                                            Grid.nfaceW(ICl,j,k));
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

          case BC_FIXED_PRESSURE :
            W[ICl-1][j][k] = W[ICl][j][k];
            W[ICl-1][j][k].p = WoW[j][k].p;
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = W[ICl][j][k] ;
            W[ICl-2][j][k].p = WoW[j][k].p;
            U[ICl-2][j][k] =  W[ICl-2][j][k].U();
            break;

          case BC_CHANNEL_INFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            //for turbulent channel flow
            // p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[ICl-1][j][k].Xc - Grid.Cell[ICl][j][k].Xc; 
            W[ICl-1][j][k] = WoW[j][k];
            W[ICl-1][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-1][j][k].p = WoW[j][k].p - dpdx*dX.x ;
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            
            dX = Grid.Cell[ICl-2][j][k].Xc - Grid.Cell[ICl][j][k].Xc;
            W[ICl-2][j][k] = WoW[j][k];
            W[ICl-2][j][k].v.x = W[ICl][j][k].v.x;
            W[ICl-2][j][k].p = WoW[j][k].p - dpdx*dX.x;
            //W[ICl-2][j][k].p = (WoW[j][k].p - 2.0/3.0*W[ICl-2][j][k].rho*W[ICl-2][j][k].k) - dpdx*dX.x;
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;
            
          case BC_PERIODIC :  
            W[ICl-1][j][k] = W[ICu-1][j][k];
            U[ICl-1][j][k] = U[ICu-1][j][k];
            W[ICl-2][j][k] = W[ICu-2][j][k];
            U[ICl-2][j][k] = U[ICu-2][j][k];
            break;

          case BC_NO_SLIP :
            W[ICl-1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[ICl][j][k],WoW[j][k], 
                                                                           Grid.nfaceW(ICl,j,k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[ICl+1][j][k], WoW[j][k], 
                                                                           Grid.nfaceW(ICl,j,k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

          case BC_MOVING_WALL :
            W[ICl-1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[ICl][j][k],WoW[j][k],
                                                                               Grid.nfaceW(ICl,j,k),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[ICl+1][j][k], WoW[j][k], 
                                                                               Grid.nfaceW(ICl,j,k),
                                                                               IPs.Moving_Wall_Velocity, 
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;
                    
          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[ICl-1][j][k] = W[ICl][j][k];
            U[ICl-1][j][k] = W[ICl-1][j][k].U();
            W[ICl-2][j][k] = W[ICl][j][k] ;
            U[ICl-2][j][k] = W[ICl-2][j][k].U();
            break;

         } /* endswitch */
         
         // Prescribe East boundary conditions.
         switch(Grid.BCtypeE[j][k]) {
          case BC_NONE :
            break;

          case BC_REFLECTION :
            W[ICu+1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[ICu][j][k],
                                                                            Grid.nfaceE(ICu,j,k));
            U[ICu+1][j][k] = W[ ICu+1][j][k].U();
            W[ICu+2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[ICu-1][j][k],
                                                                            Grid.nfaceE(ICu,j,k));
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

          case BC_FIXED_PRESSURE :
            W[ICu+1][j][k] = W[ ICu][j][k];
            W[ICu+1][j][k].p = WoE[j][k].p;
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = W[ICu-1][j][k];
            W[ICu+2][j][k].p = WoE[j][k].p; 
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

          case BC_CHANNEL_OUTFLOW:
            dpdx = IPs.Pressure_Gradient.x; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[ICu+1][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+1][j][k] = W[ICu][j][k]; 
            W[ICu+1][j][k].p = WoE[j][k].p-dpdx*dX.x;
            //W[ICu+1][j][k].p = (WoE[j][k].p-2.0/3.0*W[ICu+1][j][k].rho*W[ICu+1][j][k].k)-dpdx*dX.x;
            U[ICu+1][j][k] = W[ICu+1][j][k].U( );
            dX = Grid.Cell[ICu+2][j][k].Xc - Grid.Cell[ICu][j][k].Xc; 
            W[ICu+2][j][k] = W[ICu][j][k];
            W[ICu+2][j][k].p = WoE[j][k].p -dpdx*dX.x; 	
            //W[ICu+2][j][k].p = (WoE[j][k].p -2.0/3.0*W[ICu+2][j][k].rho*W[ICu+2][j][k].k)-dpdx*dX.x; 	
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;
            
          case BC_PERIODIC :
            W[ICu+1][j][k] = W[ICl+1][j][k];
            U[ICu+1][j][k] = U[ICl+1][j][k];
            W[ICu+2][j][k] = W[ICl+2][j][k];
            U[ICu+2][j][k] = U[ICl+2][j][k];
            break;

          case BC_NO_SLIP :
            W[ICu+1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[ICu][j][k], WoE[j][k], 
                                                                           Grid.nfaceE(ICu,j,k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[ICu-1][j][k], WoE[j][k], 
                                                                           Grid.nfaceE(ICu,j,k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;
            
          case BC_MOVING_WALL :
            W[ICu+1][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[ICu][j][k], WoE[j][k], 
                                                                               Grid.nfaceE(ICu,j,k),
                                                                               IPs.Moving_Wall_Velocity,    
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[ICu-1][j][k], WoE[j][k], 
                                                                               Grid.nfaceE(ICu,j,k),
                                                                               IPs.Moving_Wall_Velocity,   
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;
            
          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[ICu+1][j][k] = W[ICu][j][k];
            U[ICu+1][j][k] = W[ICu+1][j][k].U();
            W[ICu+2][j][k] = W[ICu][j][k];
            U[ICu+2][j][k] = W[ICu+2][j][k].U();
            break;

         } /* endswitch */
         
      } /* endfor */
   } /* endfor */
   
   for (int k = KCl - Nghost; k <= KCu + Nghost; ++k ) {
      for (int i = ICl - Nghost; i <= ICu + Nghost; ++i ) {
        // Prescribe North boundary conditions.
         switch(Grid.BCtypeN[i][k]) {
          case BC_NONE :
            break;
            
          case BC_REFLECTION :
            W[i][ JCu+1][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect( W[i][ JCu][k],
                                       Grid.nfaceN(i,  JCu,k));
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
            W[i][ JCu+2][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect( W[i][ JCu-1][k],
                                       Grid.nfaceN(i,  JCu, k));
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();
            break;
            
          case BC_FIXED_PRESSURE :
            W[i][ JCu+1][k] = W[i][ JCu][k];
            W[i][ JCu+1][k].p = WoN[i][k].p;
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
          
            W[i][ JCu+2][k] =  W[i][ JCu-1][k];
            W[i][ JCu+2][k].p = WoN[i][k].p;
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();
            break;
            
          case BC_CHANNEL_INFLOW:
            dpdy = IPs.Pressure_Gradient.y; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[i][JCu+1][k].Xc - Grid.Cell[i][JCu][k].Xc; 
            W[i][JCu+1][k].v.y = W[i][JCu][k].v.y; 
            W[i][JCu+1][k].p = WoN[i][k].p;
            //W[i][JCu+1][k].p = (WoN[i][k].p-2.0/3.0*W[i][JCu+1][k].rho*W[i][JCu+1][k].k)-dpdy*dX.y;
            
            U[i][JCu+1][k] = W[i][JCu+1][k].U( );
            
            dX = Grid.Cell[i][JCu+2][k].Xc - Grid.Cell[i][JCu][k].Xc; 
            W[i][JCu+2][k] = W[i][JCu][k];
             W[i][JCu+2][k].p = WoN[i][k].p;
            //W[i][JCu+2][k].p = (WoN[i][k].p -2.0/3.0*W[i][JCu+2][k].rho*W[i][JCu+2][k].k)-dpdy*dX.y; 	
            U[i][JCu+2][k] = W[i][JCu+2][k].U( );
            break;
            
          case BC_PERIODIC :
            W[i][JCu+1][k] = W[i][JCl+1][k];
            U[i][JCu+1][k] = U[i][JCl+1][k];
            W[i][JCu+2][k] = W[i][JCl+2][k];
            U[i][JCu+2][k] = U[i][JCl+2][k];
            break;

          case BC_NO_SLIP :
            W[i][JCu+1][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][JCu][k], WoN[i][k],
                                                                            Grid.nfaceN(i, JCu, k),
                                                                            IPs.Pressure_Gradient,
                                                                            FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U();
            W[i][JCu+2][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][JCu-1][k], WoN[i][k],
                                                                           Grid.nfaceN(i, JCu, k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;
            
          case BC_MOVING_WALL :
            W[i][JCu+1][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][JCu][k], WoN[i][k],
                                                                               Grid.nfaceN(i, JCu, k),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][JCu+1][k] = W[i][JCu+1][k].U();
            W[i][JCu+2][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][JCu-1][k], WoN[i][k],
                                                                               Grid.nfaceN(i, JCu, k),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][JCu+2][k] = W[i][JCu+2][k].U();
            break;

          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[i][ JCu+1][k] = W[i][ JCu][k];
            U[i][ JCu+1][k] =  W[i][ JCu+1][k].U();
            W[i][ JCu+2][k] =  W[i][ JCu][k];
            U[i][ JCu+2][k] =  W[i][ JCu+2][k].U();
            break;

         } /* endswitch */
    
         // Prescribe South boundary conditions.
         switch(Grid.BCtypeS[i][k]) {
          case BC_NONE :
            break;

          case BC_REFLECTION :
            W[i][JCl-1][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][JCl][k],
                                                                            Grid.nfaceS(i, JCl, k));
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            W[i][JCl-2][k] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][JCl+1][k],
                                                                            Grid.nfaceS(i, JCl, k));
            U[i][JCl-2][k] = W[i][JCl-2][k].U();
            break;
         
          case BC_FIXED_PRESSURE :
            W[i][JCl-1][k] = W[i][ JCl][k];
            W[i][JCl-1][k].p = WoS[i][k].p;
            U[i][JCl-1][k] = W[i][ JCl-1][k].U();
            W[i][JCl-2][k] = W[i][ JCl+1][k];
            W[i][JCl-2][k].p = WoS[i][k].p;
            U[i][JCl-2][k] = W[i][ JCl-2][k].U();
	    break;

          case BC_PERIODIC :
            W[i][JCl-1][k] = W[i][JCu-1][k];
            U[i][JCl-1][k] = U[i][JCu-1][k];
            W[i][JCl-2][k] = W[i][JCu-2][k];
            U[i][JCl-2][k] = U[i][JCu-2][k];
            break;

          case BC_CHANNEL_OUTFLOW:
            dpdy = IPs.Pressure_Gradient.y; 
            //for turbulent channel flow
            // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[i][JCl-1][k].Xc - Grid.Cell[i][JCl][k].Xc; 
            W[i][JCl-1][k] = WoS[i][k];
            W[i][JCl-1][k].p =  WoS[i][k].p -dpdy*dX.y;
            //W[i][JCl-1][k].p = (WoS[i][k].p )-dpdy*dX.y;
            U[i][JCl-1][k] = W[i][JCl-1][k].U();
            
            dX = Grid.Cell[i][JCl-2][k].Xc - Grid.Cell[i][JCl][k].Xc;
            W[i][JCl-2][k].p = WoS[i][k].p - dpdy*dX.y;
            //W[i][JCl-2][k].p = (WoS[i][k].p - 2.0/3.0*W[i][JCl-2][k].rho*W[i][JCl-2][k].k) - dpdy*dX.y;
            U[i][JCl-2][k] = W[i][JCl-2][k].U();
            break;

          case BC_NO_SLIP :
            W[i][JCl-1][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][JCl][k], WoS[i][k], 
                                                                           Grid.nfaceS(i, JCl,k),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] =  W[i][JCl-1][k].U();
            W[i][JCl-2][k] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][JCl+1][k], WoS[i][k],
                                                                           Grid.nfaceS(i, JCl,k),   
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] =  W[i][JCl-2][k].U();
            break;

          case BC_MOVING_WALL :
            W[i][JCl-1][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][JCl][k], WoS[i][k], 
                                                                               Grid.nfaceS(i, JCl,k),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient, 
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][JCl-1][k] =  W[i][JCl-1][k].U();
            W[i][JCl-2][k] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][JCl+1][k], WoS[i][k],
                                                                               Grid.nfaceS(i, JCl,k),
                                                                               IPs.Moving_Wall_Velocity, 
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][JCl-2][k] =  W[i][JCl-2][k].U();
            break;
            
          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[i][JCl-1][k] = W[i][ JCl][k];
            U[i][JCl-1][k] = W[i][ JCl-1][k].U();
            W[i][JCl-2][k] = W[i][ JCl][k];
            U[i][JCl-2][k] = W[i][ JCl-2][k].U();
	    break;

         } /* endswitch */

      } /* endfor */
   } /* endfor */
 
   for (int j = JCl - Nghost; j <= JCu + Nghost; ++j) {
      for (int i = ICl - Nghost; i <= ICu + Nghost; ++i) {
         // Prescribe Bottom boundary conditions.
         switch( Grid.BCtypeB[i][j]) {
          case BC_NONE :
            break;

          case BC_REFLECTION :
            W[i][j][ KCl-1] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][j][ KCl],
                                                                             Grid.nfaceBot(i, j, KCl));
            U[i][j][ KCl-1] = W[i][j][KCl-1].U();
            W[i][j][ KCl-2] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][j][ KCl+1],
                                                                             Grid.nfaceBot(i,j, KCl));
            U[i][j][ KCl-2] = W[i][j][KCl-2].U();
            break;
            
          case BC_FIXED_PRESSURE :
            W[i][j][ KCl-1] = W[i][j][ KCl];
            W[i][j][ KCl-1].p = WoB[i][j].p;
            U[i][j][ KCl-1] = W[i][j][ KCl-1].U();
            W[i][j][ KCl-2] = W[i][j][ KCl+1];
            W[i][j][ KCl-2].p = WoB[i][j].p;
            U[i][j][ KCl-2] = W[i][j][ KCl-2].U();
            break;

          case BC_CHANNEL_INFLOW:
            dpdz = IPs.Pressure_Gradient.z; 
            //for turbulent channel flow
            // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
            dX = Grid.Cell[i][j][KCl-1].Xc - Grid.Cell[i][j][KCl].Xc; 
            if ((j==JCl-1 || j==JCl-2 || j== JCu+1 || j== JCu+2 ) && 
                ((Grid.BCtypeN[i][KCl-1] == BC_NO_SLIP) ||
                 (Grid.BCtypeN[i][KCl-2] == BC_NO_SLIP) ||
                 (Grid.BCtypeS[i][KCl-1] == BC_NO_SLIP) ||
                 (Grid.BCtypeS[i][KCl-2] == BC_NO_SLIP))){
                //do not overwrite the solid no-slip wall boundary condition 
            } else {
               W[i][j][KCl-1] = WoB[i][j];
               W[i][j][KCl-1].v.z = W[i][j][KCl].v.z;
               W[i][j][KCl-1].p = WoB[i][j].p - dpdz*dX.z;
               // W[i][j][KCl-1].p = (WoB[i][j].p  - 2.0/3.0*W[i][j][KCl-1].rho* W[i][j][KCl-1].k)-dpdz*dX.z;
               U[i][j][KCl-1] = W[i][j][KCl-1].U();
               
               dX = Grid.Cell[i][j][KCl-2].Xc - Grid.Cell[i][j][KCl].Xc;
               W[i][j][KCl-2] = WoB[i][j];
               W[i][j][KCl-2].v.z = W[i][j][KCl].v.z;
               W[i][j][KCl-2].p = WoB[i][j].p - dpdz*dX.z;
               // W[i][j][KCl-2].p = (WoB[i][j].p - 2.0/3.0*W[i][j][KCl-2].rho* W[i][j][KCl-2].k ) - dpdz*dX.z;
               U[i][j][KCl-2] = W[i][j][KCl-2].U();
            } /* endif */
            break;

          case BC_PERIODIC :
            W[i][j][KCl-1] = W[i][j][KCu-1];
            U[i][j][KCl-1] = U[i][j][KCu-1];
            W[i][j][KCl-2] = W[i][j][KCu-2];
            U[i][j][KCl-2] = U[i][j][KCu-2];
            break;

          case BC_NO_SLIP :
            W[i][j][KCl-1] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][j][KCl], WoB[i][j],
                                                                           Grid.nfaceBot(i, j, KCl),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] = W[i][j][KCl-1].U(); 
            W[i][j][KCl-2] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][j][KCl+1], WoB[i][j],
                                                                           Grid.nfaceBot(i, j, KCl),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;
            
          case BC_MOVING_WALL :
            W[i][j][KCl-1] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][j][KCl], WoB[i][j],
                                                                               Grid.nfaceBot(i, j,KCl),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-1] = W[i][j][KCl-1].U();
            W[i][j][KCl-2] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][j][KCl+1], WoB[i][j],
                                                                               Grid.nfaceBot(i, j ,KCl),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][j][KCl-2] = W[i][j][KCl-2].U();
            break;

          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[i][j][ KCl-1] = W[i][j][ KCl];
            U[i][j][ KCl-1] = W[i][j][ KCl-1].U();
            W[i][j][ KCl-2] = W[i][j][ KCl];
            U[i][j][ KCl-2] = W[i][j][ KCl-2].U();
            break;
            
         } /* endswitch */
         
         // Prescribe Top boundary conditions.
         switch(Grid.BCtypeT[i][j]) {
          case BC_NONE :
            break;
            
          case BC_REFLECTION :
            W[i][j][KCu+1] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][j][ KCu],
                                                                            Grid.nfaceTop(i, j, KCu));
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = FANS3D_ThermallyPerfect_KOmega_pState::Reflect(W[i][j][ KCu-1],
                                                                            Grid.nfaceTop(i, j, KCu));
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;
          case BC_FIXED_PRESSURE :
            W[i][j][KCu+1] = W[i][j][ KCu];
            W[i][j][KCu+1].p = WoT[i][j].p;
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = W[i][j][KCu-1];
            W[i][j][KCu+2].p = WoT[i][j].p;
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

          case BC_CHANNEL_OUTFLOW:
            dpdz = IPs.Pressure_Gradient.z; 
            // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
            dX = Grid.Cell[i][j][KCu+1].Xc - Grid.Cell[i][j][KCu].Xc; 
            W[i][j][KCu+1] = W[i][j][KCu]; 
            W[i][j][KCu+1].p = (WoT[i][j].p)-dpdz*dX.z;
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            dX = Grid.Cell[i][j][KCu+2].Xc - Grid.Cell[i][j][KCu].Xc; 
            W[i][j][KCu+2] = W[i][j][KCu];
            W[i][j][KCu+2].p = (WoT[i][j].p)-dpdz*dX.z; 	
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

          case BC_PERIODIC :
            W[i][j][KCu+1] = W[i][j][KCl+1];
            U[i][j][KCu+1] = U[i][j][KCl+1];
            W[i][j][KCu+2] = W[i][j][KCl+2];
            U[i][j][KCu+2] = U[i][j][KCl+2];
            break;

          case BC_NO_SLIP :
            W[i][j][KCu+1] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][j][KCu], WoT[i][j],
                                                                           Grid.nfaceTop(i, j , KCu),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(W[i][j][KCu -1], WoT[i][j],
                                                                           Grid.nfaceTop(i, j , KCu),
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;
            
          case BC_MOVING_WALL :
            W[i][j][KCu+1] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][j][KCu], WoT[i][j],
                                                                               Grid.nfaceTop(i, j , KCu),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(W[i][j][KCu -1], WoT[i][j],
                                                                               Grid.nfaceTop(i, j , KCu),
                                                                               IPs.Moving_Wall_Velocity,
                                                                               IPs.Pressure_Gradient,
                                                                               FIXED_TEMPERATURE_WALL);
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;
 
          case BC_CONSTANT_EXTRAPOLATION :
	  default :
            W[i][j][KCu+1] = W[i][j][KCu];
            U[i][j][KCu+1] = W[i][j][KCu+1].U();
            W[i][j][KCu+2] = W[i][j][KCu];
            U[i][j][KCu+2] = W[i][j][KCu+2].U();
            break;

         } /* endswitch */

      } /* endfor */
   } /* endfor */
               
   // compute y+ etc.;
   Wall_Shear();
   
}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for the    *
 * specified hexa solution block according to           *
 * the Courant-Friedrichs-Lewy condition.               *
 *                                                      *
 ********************************************************/
template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                  FANS3D_ThermallyPerfect_KOmega_cState>::
CFL(Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs){
   
   double dtMin, d_i, d_j, d_k, v_i, v_j, v_k, a, dt_vis, nv, nv_t;
   
   dtMin = MILLION;
   
   for (int k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) {
      for (int j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         for (int i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = ZERO;
            } else {
               d_i = TWO*(Grid.volume(i,j,k)/
                          (Grid.AfaceE(i, j, k)+ Grid.AfaceW(i, j, k)));
               d_j = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceN(i, j, k)+ Grid.AfaceS(i, j, k)));
               d_k = TWO*( Grid.volume(i, j, k)/
                           (Grid.AfaceTop(i, j, k)+ Grid.AfaceBot(i, j, k)));
               v_i = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           (Grid.nfaceE(i, j, k)- Grid.nfaceW(i, j, k)));
               v_j = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                           ( Grid.nfaceN(i, j, k)- Grid.nfaceS(i, j, k)));
               v_k = HALF*(W[i][j][k].rhov()/W[i][j][k].rho*
                            (Grid.nfaceTop(i, j, k)- Grid.nfaceBot(i, j, k)));
               a =  W[i][j][k].a();
               dt[i][j][k] = min(min(d_i/(a+fabs(v_i)), 
                                     d_j/(a+fabs(v_j))),
                                     d_k/(a+fabs(v_k)));
               
               if (Flow_Type != FLOWTYPE_INVISCID) {  
                  nv = W[i][j][k].mu()/W[i][j][k].rho;
         
                  if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
                      Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
                     nv_t = W[i][j][k].mu_t()/W[i][j][k].rho; 
                     nv = max(nv, nv_t);
                  } /* endif */
        
                  dt_vis = min(min((d_i*d_i)/(3.0*nv), 
                                   (d_j*d_j)/(3.0*nv)), 
                                   (d_k*d_k)/(3.0*nv)); 
                  dt[i][j][k]  = min(dt_vis, dt[i][j][k]);
               } /* endif */
               
               dtMin = min(dtMin,  dt[i][j][k]);
               
            } /* endif */
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   for (int k  =  KCl- Nghost ; k <=  KCu+ Nghost ; ++k ) {
      for (int j  =  JCl- Nghost ; j <=  JCu+ Nghost ; ++j ) {
         for (int i =  ICl- Nghost ; i <=  ICu+ Nghost ; ++i ) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = dtMin;
            } /* endif */
         } /* endfor */
      } /* endfor */
   } /* endfor */
   
   /* Return the global time step. */
   
   return (dtMin);

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine determines the solution residuals for a *
 * given stage of a variety of multi-stage explicit     *
 * time integration schemes for a given solution block. *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                                          FANS3D_ThermallyPerfect_KOmega_cState> &IPs) {

   int i, j, k,  k_residual;
   double omega;
   Vector3D dX;
   
   FANS3D_ThermallyPerfect_KOmega_pState Wl, Wr;
   FANS3D_ThermallyPerfect_KOmega_cState Flux, Temp;

   FANS3D_ThermallyPerfect_KOmega_cState U_VACUUM;
   U_VACUUM.Vacuum();
   FANS3D_ThermallyPerfect_KOmega_pState W_VACUUM;
   W_VACUUM.Vacuum();
   Wl = W_VACUUM;
   Wr = W_VACUUM;
   Flux = U_VACUUM;
   Temp = U_VACUUM;


  /* Evaluate the time step fraction and residual storage location 
   *     *       for the stage. */
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
         } /* endif */
      } /* endif */
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
   
   
   /* Perform the linear reconstruction within each cell
      of the computational grid for this stage. */
   
   switch(IPs.i_Reconstruction) {
      
   case RECONSTRUCTION_LEAST_SQUARES :

      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);

      break;
   default:
      Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
      break;
   } /* endswitch */

   // compute y+ etc.;
   Wall_Shear( );
   
   /* Evaluate the time rate of change of the solution
      (i.e., the solution residuals) using a second-order
      limited upwind scheme with a variety of flux functions. */
   
   // Add i-direction (zeta-direction) fluxes.
   for ( k  =  KCl-1 ; k <=  KCu+1 ; ++k ){
      
      for ( j  =  JCl-1 ; j <=  JCu+1 ; ++j ) {
         if ( i_stage == 1 ) {
            Uo[ ICl-1][j][k] =  U[ ICl-1][j][k];
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
         } else {
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
         } /* endif */
         
         for ( i =  ICl-1 ; i <=  ICu ; ++i ) {
            if ( i_stage == 1 ) {
               Uo[i+1][j][k] =  U[i+1][j][k];
               dUdt[i+1][j][k][k_residual] = U_VACUUM;
            } else if (( j >  JCl-1 && j <  JCu+1 ) 
                       && (k >  KCl-1 && k <  KCu+1) ){
               switch(IPs.i_Time_Integration) {
               case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
                  //dUdt[i+1][j][k][k_residual] = 
                  //  dUdt[i+1][j][k][k_residual];
                  break;
               case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
		  if (IPs.N_Stage == 2) {
                     //dUdt[i+1][j][k_residual] = 
                     // dUdt[i+1][j][k_residual];
                  } else if (IPs.N_Stage == 4 && i_stage == 4) {
                     dUdt[i+1][j][k][k_residual] = 
                        dUdt[i+1][j][k][0] + 
                        TWO*dUdt[i+1][j][k][1] +TWO*dUdt[i+1][j][k][2];
                  } else {
         	     dUdt[i+1][j][k][k_residual].Vacuum();
                  } /* endif */
                  break;

               case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  break;
               default:
                  dUdt[i+1][j][k][k_residual] = U_VACUUM;
                  break;
               } /* endswitch */
            } /* endif */
            if (( j >  JCl-1 && j <  JCu+1 ) && 
                ( k >  KCl-1 && k <  KCu+1 )) {
               
               /* Evaluate the cell interface i-direction fluxes. */
               if (i ==  ICl-1 &&  (Grid.BCtypeW[j][k] == BC_REFLECTION ||
                                    Grid.BCtypeW[j][k] == BC_NO_SLIP||
                                    Grid.BCtypeW[j][k] == BC_MOVING_WALL)  ) {
                               
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  
                  Wr =  W[i+1][j][k] + 
                     ( phi[i+1][j][k]^dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^dWdz[i+1][j][k])*dX.z;
                                 
                  if ( Grid.BCtypeW[j][k] == BC_REFLECTION) {
                     Wl =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wr,  Grid.nfaceW(i+1, j, k));
                  }
                  if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                     Wl =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k), 
                                                                          IPs.Pressure_Gradient,
                                                                          FIXED_TEMPERATURE_WALL);
                  }
                  if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                     Wl =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wr, WoW[j][k], Grid.nfaceW(i+1, j, k),
                                                                              IPs.Moving_Wall_Velocity,
                                                                              IPs.Pressure_Gradient,
                                                                              FIXED_TEMPERATURE_WALL);
                  }
               } else if (i ==  ICu && 
                          ( Grid.BCtypeE[j][k] == BC_REFLECTION||
                            Grid.BCtypeE[j][k] == BC_NO_SLIP||
                            Grid.BCtypeE[j][k] == BC_MOVING_WALL)) {
                  
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  
                  if ( Grid.BCtypeE[j][k] == BC_REFLECTION) {
                     Wr =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wl,  Grid.nfaceE(i, j, k));
                     
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                     Wr =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wl, WoE[j][k], Grid.nfaceE(i, j, k), 
                                                                          IPs.Pressure_Gradient, 
                                                                          FIXED_TEMPERATURE_WALL);
                     
                  } 
                  if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                     Wr =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wl, WoE[j][k], Grid.nfaceE(i, j, k),  
                                                                              IPs.Moving_Wall_Velocity,
                                                                              IPs.Pressure_Gradient,
                                                                              FIXED_TEMPERATURE_WALL);
                  }
                  
               } else { 
                  dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
                  Wl =  W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               
                  dX =  Grid.xfaceW(i+1, j, k)- Grid.Cell[i+1][j][k].Xc;
                  Wr =  W[i+1][j][k] + 
                     ( phi[i+1][j][k]^ dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^ dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^ dWdz[i+1][j][k])*dX.z;
               } /* endif */
               
               switch(IPs.i_Flux_Function) {
               case FLUX_FUNCTION_HLLE :
                  Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                  break;
               case FLUX_FUNCTION_ROE :
                  Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceE(i, j, k));
                  break;
               } /* endswitch */

               // add viscous flux in i direction
               Flux -=  FANS3D_ThermallyPerfect_KOmega_pState::FluxViscous_n(Wl, Wr,W[i][j][k], W[i+1][j][k], 
                                                                             dWdx[i][j][k], dWdy[i][j][k],dWdz[i][j][k],
                                                                             dWdx[i+1][j][k], dWdy[i+1][j][k],dWdz[i+1][j][k],
                                                                             Grid.nfaceE(i, j, k), Grid.Voe(i, j, k),  
                                                                             Grid.delta_oe(i, j, k),
                                                                             Grid.volume(i, j, k), 
                                                                            Grid.volume(i+1, j, k));

               /* Evaluate cell-averaged solution changes. */
               
               dUdt[i][j][k][k_residual] -= 
                  (IPs.CFL_Number* dt[i][j][k])*
                  Flux* Grid.AfaceE(i, j, k)/
                  Grid.volume(i, j, k);
               dUdt[i+1][j][k][k_residual] += 
                  (IPs.CFL_Number* dt[i+1][j][k])*
                  Flux* Grid.AfaceW(i+1, j, k)/
                  Grid.volume(i+1, j ,k);
                                 
               /* Include source terms associated with turbulence model */
               dUdt[i][j][k][k_residual] += (IPs.CFL_Number* dt[i][j][k])*
                  FANS3D_ThermallyPerfect_KOmega_pState::Sturbulence(W[i][j][k], 
                                                                     dWdx[i][j][k], 
                                                                     dWdy[i][j][k], 
                                                                     dWdz[i][j][k]);
     
	   /* Save west and east face boundary flux. */
               
               //    if (i ==  ICl-1) {
//                     FluxW[j] = -Flux* Grid.lfaceW(i+1, j);
//                 } else if (i ==  ICu) {
//                     FluxE[j] = Flux* Grid.lfaceE(i, j);
//                 } /* endif */ 
               
            } /* endif */
         } /* endfor */
          
         if (( j >  JCl-1 && j <  JCu+1 ) && 
             ( k >  KCl-1 && k <  KCu+1 ) ) {
            dUdt[ ICl-1][j][k][k_residual] = U_VACUUM;
            dUdt[ ICu+1][j][k][k_residual] = U_VACUUM;
         } /* endif */
      } /* endfor */
   } /* endfor */
 
   // Add j-direction (eta-direction) fluxes.
   for ( k =  KCl ; k <=  KCu ; ++k ) {
      
      for ( i =  ICl ; i <=  ICu ; ++i ) {
         for ( j  =  JCl-1 ; j <=  JCu ; ++j ) {
            
            /* Evaluate the cell interface j-direction fluxes. */
            if (j ==  JCl-1 && 
                ( Grid.BCtypeS[i][k] == BC_REFLECTION ||
                  Grid.BCtypeS[i][k] == BC_NO_SLIP||
                  Grid.BCtypeS[i][k] == BC_MOVING_WALL)) {
               
               dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
               Wr =  W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y+
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
               if ( Grid.BCtypeS[i][k] == BC_REFLECTION) {
                  Wl =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wr,  Grid.nfaceS(i, j+1, k));
                }
                if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                   Wl =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                                                        IPs.Pressure_Gradient,
                                FIXED_TEMPERATURE_WALL);
                }
                if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                   Wl =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wr, WoS[i][k], Grid.nfaceS(i, j+1, k),
                                                                            IPs.Moving_Wall_Velocity,
                                                                            IPs.Pressure_Gradient,
                                                                            FIXED_TEMPERATURE_WALL);
                }

             } else if (j ==  JCu && 
                        ( Grid.BCtypeN[i][k] == BC_REFLECTION ||
                          Grid.BCtypeN[i][k] == BC_NO_SLIP||
                          Grid.BCtypeN[i][k] == BC_MOVING_WALL)) {
                
                dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                if ( Grid.BCtypeN[i][k] == BC_REFLECTION) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wl,  Grid.nfaceN(i, j, k));
                }
                if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                                                        IPs.Pressure_Gradient,
                                                                        FIXED_TEMPERATURE_WALL );
                }
                if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wl, WoN[i][k], Grid.nfaceN(i, j, k),
                                                                            IPs.Moving_Wall_Velocity,
                                                                            IPs.Pressure_Gradient,
                                                                            FIXED_TEMPERATURE_WALL );
                }
             } else {
                dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX =  Grid.xfaceS(i, j+1, k)- Grid.Cell[i][j+1][k].Xc;
                Wr =  W[i][j+1][k] +
                   ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                   ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y +
                   ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;
          
                        
          } /* endif */

             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceN(i, j, k));
                break;
             case FLUX_FUNCTION_ROE :
                Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceN(i, j, k));
                break;
             } /* endswitch */
              
             // add viscous flux in j direction
             
             Flux -=  FANS3D_ThermallyPerfect_KOmega_pState::FluxViscous_n(Wl, Wr, W[i][j][k], W[i][j+1][k],
                                                                           dWdx[i][j][k], dWdy[i][j][k],dWdz[i][j][k],
                                                                           dWdx[i][j+1][k], dWdy[i][j+1][k],dWdz[i][j+1][k],
                                                                           Grid.nfaceN(i, j, k),  Grid.Von(i, j, k),
                                                                           Grid.delta_on(i, j, k), Grid.volume(i, j, k), 
                                                                           Grid.volume(i, j+1, k));
             
             /* Evaluate cell-averaged solution changes. */
             dUdt[i][j][k][k_residual] -= 
                (IPs.CFL_Number* dt[i][j][k])*
                Flux* Grid.AfaceN(i, j, k)/
                Grid.volume(i, j, k);
              dUdt[i][j+1][k][k_residual] += 
                (IPs.CFL_Number* dt[i][j+1][k])*
                Flux* Grid.AfaceS(i, j+1, k)/
                 Grid.volume(i, j+1, k);
    
              
          /* Save south and north face boundary flux. */

          // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */

          } /* endfor */
           dUdt[i][ JCl-1][k][k_residual] = U_VACUUM;
           dUdt[i][ JCu+1][k][k_residual] = U_VACUUM;
       } /* endfor */
   } /* endfor */

   // Add k-direction (gamma-direction) fluxes.
   for ( i =  ICl; i <=  ICu ; ++i ) {
   
      for ( j  =  JCl ; j <=  JCu ; ++j ){
         for ( k =  KCl-1 ; k <=  KCu ; ++k )  {
            
            /* Evaluate the cell interface j-direction fluxes. */
            if (k ==  KCl-1 && 
                ( Grid.BCtypeB[i][j] == BC_REFLECTION  || 
                  Grid.BCtypeB[i][j] == BC_NO_SLIP ||
                  Grid.BCtypeB[i][j] == BC_MOVING_WALL)) {
                
               dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr =  W[i][j][k+1] +
                  ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                  ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y+
                  ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
               if ( Grid.BCtypeB[i][j] == BC_REFLECTION) {
                  Wl =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wr,  Grid.nfaceBot(i, j, k+1));
               }
               if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                  Wl =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1), 
                                                                       IPs.Pressure_Gradient,
                                                                       FIXED_TEMPERATURE_WALL);
               }
               if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                  Wl =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wr, WoB[i][j], Grid.nfaceBot(i, j, k+1),
                                                                           IPs.Moving_Wall_Velocity,
                                                                           IPs.Pressure_Gradient,
                                                                           FIXED_TEMPERATURE_WALL);
               }


            } else if (k ==  KCu && 
                       ( Grid.BCtypeT[i][j] == BC_REFLECTION ||
                         Grid.BCtypeT[i][j] == BC_NO_SLIP ||
                         Grid.BCtypeT[i][j] == BC_MOVING_WALL)) {
                
               dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                if ( Grid.BCtypeT[i][j] == BC_REFLECTION) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::Reflect(Wl,  Grid.nfaceTop(i, j, k));
                }
                if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::NoSlip(Wl, WoT[i][j], Grid.nfaceTop(i, j, k), 
                                                                       IPs.Pressure_Gradient,
                                                                        FIXED_TEMPERATURE_WALL );
                }
                if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                   Wr =  FANS3D_ThermallyPerfect_KOmega_pState::MovingWall(Wl, WoT[i][j], Grid.nfaceTop(i, j, k),
                                                                            IPs.Moving_Wall_Velocity,
                                                                            IPs.Pressure_Gradient, 
                                                                            FIXED_TEMPERATURE_WALL );
                }

             } else {
                dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
                Wl =  W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                dX =  Grid.xfaceBot(i, j, k+1)- Grid.Cell[i][j][k+1].Xc;
                Wr =  W[i][j][k+1] +
                   ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                   ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y +
                   ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;
             } /* endif */
             
             switch(IPs.i_Flux_Function) {
             case FLUX_FUNCTION_HLLE :
                Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
                break;
             case FLUX_FUNCTION_ROE :
                Flux =  FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_n(Wl, Wr,  Grid.nfaceTop(i, j, k));
                break;
             } /* endswitch */
             
             // add viscous flux in k direction
             Flux -=  FANS3D_ThermallyPerfect_KOmega_pState::FluxViscous_n(Wl, Wr,W[i][j][k], W[i][j][k+1], 
                                   dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                   dWdx[i][j][k+1], dWdy[i][j][k+1], dWdz[i][j][k+1],
                                   Grid.nfaceTop(i, j, k), Grid.Vot(i, j, k), 
                                   Grid.delta_ot(i, j, k), Grid.volume(i, j, k), Grid.volume(i, j, k+1));

             /* Evaluate cell-averaged solution changes. */
             
             dUdt[i][j][k][k_residual] -= 
                (IPs.CFL_Number* dt[i][j][k])*
                Flux* Grid.AfaceTop(i, j, k)/
                Grid.volume(i, j, k);
             dUdt[i][j][k+1][k_residual] += 
                (IPs.CFL_Number* dt[i][j][k+1])*
                Flux* Grid.AfaceBot(i, j, k+1)/
                Grid.volume(i, j, k+1);
             
             /* Save top and bottom face boundary flux. */
             
             // if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
             
          } /* endfor */
          dUdt[i][j][ KCl-1][k_residual] = U_VACUUM;
          dUdt[i][j][ KCu+1][k_residual] = U_VACUUM;

      } /* endfor */

   } /* endfor */
   
   // Residual "correction" for sublayer and log-layer based on 
   // wall boundary treatment methods.
   if (IPs.Turbulence_IP.Near_Wall_Boundary_Treatment == 0){
      for ( k  =  KCl ; k <=  KCu ; ++k ) 
         for ( j  =  JCl ; j <=  JCu ; ++j ) 
            for ( i =  ICl ; i <=  ICu ; ++i ) {
               U[i][j][k].k_omega_model.automatic_wall_treatment_residual_evaluaton(
                  i, j, k, WallData[i][j][k], Grid,  dUdt[i][j][k][k_residual].rhok, dUdt[i][j][k][k_residual].rhoomega);
            }
      
   }else if (IPs.Turbulence_IP.Near_Wall_Boundary_Treatment == 2){
      for ( k  =  KCl ; k <=  KCu ; ++k ) 
         for ( j  =  JCl ; j <=  JCu ; ++j ) 
            for ( i =  ICl ; i <=  ICu ; ++i ) {
               U[i][j][k].k_omega_model.low_Reynolds_number_formulation_residual_evaluaton(
                  i, j, k, WallData[i][j][k], Grid, dUdt[i][j][k][k_residual].rhoomega);
            }
      
   }else{
      for ( k  =  KCl ; k <=  KCu ; ++k ) 
         for ( j  =  JCl ; j <=  JCu ; ++j ) 
            for ( i =  ICl ; i <=  ICu ; ++i ) {
               U[i][j][k].k_omega_model.wall_function_residual_evaluation(
                  WallData[i][j][k], dUdt[i][j][k][k_residual].rhok,
                  dUdt[i][j][k][k_residual].rhoomega);
            }
      
   } // end of the residual "correction".
   
   
   /* Residual for the stage successfully calculated. */
   
   
   return (0);
    
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
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
               FANS3D_ThermallyPerfect_KOmega_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState,
                                                     FANS3D_ThermallyPerfect_KOmega_cState> &IPs){
   
   int k_residual;
   double omega;
   
   int num_vars = NumVar();

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
         } /* endif */
      } /* endif */
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

   for (int k =  KCl ; k <= KCu ; ++k ) {
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
            
            // Wall function and low Reynolds number formulation  
            if (IPs.Turbulence_IP.Near_Wall_Boundary_Treatment == 0) {
               // automatic wall treatment formulation
               U[i][j][k].k_omega_model.automatic_wall_treatment(U[i][j][k].rho, 
                                                                 W[i][j][k].mu(), 
                                                                 i, j, k, 
                                                                 WallData[i][j][k],
                                                                 Grid,
                                                                 U[i][j][k].rhok, 
                                                                 U[i][j][k].rhoomega);
                       
            } else if(IPs.Turbulence_IP.Near_Wall_Boundary_Treatment == 2) {
               // low Reynolds number formulation
               U[i][j][k].k_omega_model.low_Reynolds_number_formulation(
                  U[i][j][k].rho, W[i][j][k].mu(), 
                  i, j, k, WallData[i][j][k],
                  Grid, U[i][j][k].rhoomega);
               
            } else { // wall function
               U[i][j][k].k_omega_model.wall_function(U[i][j][k].rho, 
                                                      W[i][j][k].mu(), 
                                                      WallData[i][j][k], 
                                                      U[i][j][k].rhok, 
                                                      U[i][j][k].rhoomega); 
            } /* endif */
            
            // Check physical validity of update solution state
            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Mass Fractions, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] 
                     << " omega = " << omega << "\n";
		return (1);
              } /* endif */

            } else {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Mass Fractions, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] 
                     << " omega = " << omega << "\n";
		return (1);
              } /* endif */

            } /* endif */
            
            W[i][j][k] = U[i][j][k].W();
            
         } /* endfor */    	 
      } /* endfor */    
   } /* endfor */
    
   /* Solution successfully updated. */
    
   return (0);   
      
}

/**********************************************************************
 * Routine: Wall_Shear                                                *
 *                                                                    *
 * This routine determines the wall shear stress, friction velocity,  *
 * and yplus at each cell given the normal distance to the wall.      *
 *                                                                    *
 **********************************************************************/
template<> 
int Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState,  
               FANS3D_ThermallyPerfect_KOmega_cState>::Wall_Shear(void) { 
   
   for (int k = Grid.KCl- Grid.Nghost ; k <= Grid.KCu+ Grid.Nghost ; ++k )
      for (int j = Grid.JCl- Grid.Nghost ; j <= Grid.JCu+ Grid.Nghost ; ++j ) 
         for (int i = Grid.ICl- Grid.Nghost ; i <= Grid.ICu+ Grid.Nghost ; ++i ) {
            
            // for ghost cells close to solid wall, the ywall is set to be zero,
            // this is related to the application of wall function (see the logic in the wall fucntion code
            // it may be not necessary, since the noslip boundary condition, tangiential velocity
            // is zero, so the wall function will most likely to be applied correctly.
            // however, it is not over-excessive to have these boundary checks.
            
            // Check West boundary.
            if (Grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
                Grid.BCtypeW[j][k] == BC_NO_SLIP  ||
                Grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeW[j][k] == BC_ADIABATIC_WALL ||
                Grid.BCtypeW[j][k] == BC_BURNING_SURFACE){
               
               WallData[ Grid.ICl-2][j][k].ywall = ZERO;
               WallData[ Grid.ICl-1][j][k].ywall = ZERO;
            } /* endif */
            
            // Check East boundary.
            if (Grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
                Grid.BCtypeE[j][k] == BC_NO_SLIP  ||
                Grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeE[j][k] == BC_ADIABATIC_WALL ||
                Grid.BCtypeE[j][k] == BC_BURNING_SURFACE) {
               WallData[ Grid.ICu+1][j][k].ywall = ZERO;
               WallData[ Grid.ICu+2][j][k].ywall = ZERO;
               
            } /* endif */
            
            // Check South boundary.
            if (Grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
                Grid.BCtypeS[i][k] == BC_NO_SLIP  ||
                Grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeS[i][k] == BC_ADIABATIC_WALL ||
                Grid.BCtypeS[i][k] == BC_BURNING_SURFACE) {
               
               WallData[i][ Grid.JCl-2][k].ywall = ZERO;
               WallData[i][ Grid.JCl-1][k].ywall = ZERO;
            } /* endif */
               
               
            // Check North boundary.
            if (Grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
                Grid.BCtypeN[i][k] == BC_NO_SLIP  ||
                Grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeN[i][k] == BC_ADIABATIC_WALL ||
                Grid.BCtypeN[i][k] == BC_BURNING_SURFACE) {
               
               WallData[i][ Grid.JCu+1][k].ywall = ZERO;
               WallData[i][ Grid.JCu+2][k].ywall = ZERO;
                  
            } /* endif */
               
            // Check Bottom boundary.
            
            if (Grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
                Grid.BCtypeB[i][j] == BC_NO_SLIP  ||
                Grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeB[i][j] == BC_ADIABATIC_WALL ||
                Grid.BCtypeB[i][j] == BC_BURNING_SURFACE) {
               
               WallData[i][j][ Grid.KCl-1].ywall = ZERO;
               WallData[i][j][ Grid.KCl-2].ywall = ZERO;
            }
            
            // Check Top boundary.
            if (Grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
                Grid.BCtypeT[i][j] == BC_NO_SLIP  ||
                Grid.BCtypeT[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeT[i][j] == BC_ADIABATIC_WALL ||
                Grid.BCtypeT[i][j] == BC_BURNING_SURFACE) {
               
               WallData[i][j][ Grid.KCu+1].ywall = ZERO;
               WallData[i][j][ Grid.KCu+2].ywall = ZERO;
               
            }
               
            WallData[i][j][k].utau = Wall_Friction_Velocity(i, j, k);//xinfeng: note: use the experimental data. utau_experimental = 0.3.
            WallData[i][j][k].tauw =  W[i][j][k].rho * sqr(WallData[i][j][k].utau);
            WallData[i][j][k].yplus =  WallData[i][j][k].utau*WallData[i][j][k].ywall/
               ( W[i][j][k].mu()/ W[i][j][k].rho);
            
                             
         } /* endfor */

   return 0;
      
}

/**********************************************************************
 * Routine: Friction_Velocity                                         *
 *                                                                    *
 * This routine determines the friction velocity.                     *
 *                                                                    *
 **********************************************************************/
template<>
double Hexa_Block<FANS3D_ThermallyPerfect_KOmega_pState, 
                  FANS3D_ThermallyPerfect_KOmega_cState>::
Wall_Friction_Velocity(const int i, 
                       const int j, 
                       const int k) {
 
   if( Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA){
      exit(1);
   } else {
      double y, yplus, yplus_o;
      double kappa, Const, Betastar, nu, rho_wall;
      double u_t0, epsilon, Friction_Velocity,value, f0,  
         tangentialvelocity, normalvelocity;
      double q, u_t1;
      int n, m;
      
      y =  WallData[i][j][k].ywall;
      
      kappa = W[i][j][k].k_omega_model.Karman_const;
      Const = W[i][j][k].k_omega_model.C_const;
      Const = exp(kappa*Const);
      Betastar = W[i][j][k].k_omega_model.beta_star;
      nu = W[i][j][k].mu()/W[i][j][k].rho;
      tangentialvelocity = abs(W[i][j][k].v);
           
      // tangentialvelocity = abs(W[i][j][k].v - (W[i][j][k].v * WallData[i][j][k].nwall)*WallData[i][j][k].nwall);
      // note: may 31, 2005 ; the total velocity seems work better for bluff body configuration
      if(y==ZERO){// for those ghost cells close to the solid walls
         Friction_Velocity = ZERO;
      }else{
         Friction_Velocity = sqrt(nu*tangentialvelocity/y);
     
      }
      
      yplus = y*Friction_Velocity/nu;
      yplus_o = 10.80487081;
      //solving from maple by getting the crosspoint of u+ = y+ ----
      //u+ = 1/kappa ln(y+)+C;
      
      if (yplus<yplus_o){
         //* some flow geometry type or some flow field type,
         //there are some regions that utangential velocity ends up to be ZERO*//
         // in which case, the zero utau is problematic for the wall function, so the following check is necessary.        //
         if (tangentialvelocity ==ZERO &&
             ( (Grid.BCtypeW[j][k] != BC_WALL_VISCOUS ||
                Grid.BCtypeW[j][k] != BC_NO_SLIP  ||
                Grid.BCtypeW[j][k] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeW[j][k] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeW[j][k] != BC_ADIABATIC_WALL ||
                Grid.BCtypeW[j][k] != BC_BURNING_SURFACE) ||
               (Grid.BCtypeE[j][k] != BC_WALL_VISCOUS ||
                Grid.BCtypeE[j][k] != BC_NO_SLIP  ||
                Grid.BCtypeE[j][k] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeE[j][k] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeE[j][k] != BC_ADIABATIC_WALL ||
                Grid.BCtypeE[j][k] != BC_BURNING_SURFACE) ||
               (Grid.BCtypeS[i][k] != BC_WALL_VISCOUS ||
                Grid.BCtypeS[i][k] != BC_NO_SLIP  ||
                Grid.BCtypeS[i][k] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeS[i][k] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeS[i][k] != BC_ADIABATIC_WALL ||
                Grid.BCtypeS[i][k] != BC_BURNING_SURFACE) ||
               (Grid.BCtypeN[i][k] != BC_WALL_VISCOUS ||
                Grid.BCtypeN[i][k] != BC_NO_SLIP  ||
                Grid.BCtypeN[i][k] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeN[i][k] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeN[i][k] != BC_ADIABATIC_WALL ||
                Grid.BCtypeN[i][k] != BC_BURNING_SURFACE) ||
               (Grid.BCtypeB[i][j] != BC_WALL_VISCOUS ||
                Grid.BCtypeB[i][j] != BC_NO_SLIP  ||
                Grid.BCtypeB[i][j] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeB[i][j] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeB[i][j] != BC_ADIABATIC_WALL ||
                Grid.BCtypeB[i][j] != BC_BURNING_SURFACE) ||
               (Grid.BCtypeT[i][j] != BC_WALL_VISCOUS ||
                Grid.BCtypeT[i][j] != BC_NO_SLIP  ||
                Grid.BCtypeT[i][j] != BC_WALL_VISCOUS_ISOTHERMAL ||
                Grid.BCtypeT[i][j] != BC_WALL_VISCOUS_HEATFLUX ||
                Grid.BCtypeT[i][j] != BC_ADIABATIC_WALL ||
                Grid.BCtypeT[i][j] != BC_BURNING_SURFACE))){
            
            // normalvelocity = abs((W[i][j][k].v * WallData[i][j][k].nwall)*WallData[i][j][k].nwall);
            // Friction_Velocity = (5.0/100.0*normalvelocity)*pow(W[i][j][k].beta_star, 0.25);
            Friction_Velocity = 5.0/100.0*abs(W[i][j][k].v)*pow(W[i][j][k].k_omega_model.beta_star, 0.25);
            // cout<<"\n  fix :  Friction_Velocity= "<< Friction_Velocity<<endl;
         }
         
         //cout<<"\n y+<y+o :   Friction_Velocity= "<< Friction_Velocity<<endl;
         return  (Friction_Velocity);
      
      }else{
         //iteratively solve to obtain the friction velocity
         //initial guess with an "analytical" solution
         u_t0 = 1.1;
         n = 0;
         m = 20;
         epsilon= 1e-3;
         
         nu = W[i][j][k].mu()/W[i][j][k].rho;
         value = log(Const*y*u_t0/nu);
         f0 = u_t0 - kappa*tangentialvelocity/value;
         
         do {
            value = log(Const*y*u_t0/nu);
            q = 1.0+ kappa*tangentialvelocity/(u_t0*value*value);
            u_t1 = u_t0 - f0/q;
            u_t0 = u_t1;
            value = log(Const*y*u_t0/nu);
            f0 = u_t0 - kappa*tangentialvelocity/value;
            n=n+1;
            
         } while((fabs(f0)>= epsilon) && (n<=m));
         
         Friction_Velocity = u_t1;
         
         return  Friction_Velocity ;
      }
   }//end of turbulent case 
   
} 

 

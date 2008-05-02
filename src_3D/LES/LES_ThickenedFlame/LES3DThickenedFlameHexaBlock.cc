
#ifndef  _LES3DTF_HEXA_BLOCK_INCLUDED
#include "LES3DThickenedFlameHexaBlock.h"
#endif // _LES3DTF_HEXA_BLOCK_INCLUDED


/******************************************************************
 * Hexa_Block::NumVar -- Return the number of solution variables. *
 ******************************************************************/
template<>
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::NumVar() {
  return (W[0][0][0].num_vars+2);
}

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
void Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
Output_Tecplot(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs,
               const int Number_of_Time_Steps,
               const double &Time,
               const int Block_Number,
               const int Output_Title,
               ostream &Out_File) {
 
   LES3DTF_pState W_node;

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);

   
   /* Output node solution data. */
   
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
               << "\"k\" \\ \n";
      //n species mass fractions names
      for (int i =0 ; i < W[0][0][0].ns ; ++i){
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      } /* endfor */
      
      Out_File <<"\"WF\" \\ \n" 
	       <<"\"TF\" \\ \n"
	       <<"\"Iso_c=0.5\" \\ \n";

      Out_File <<"\"T\" \\ \n"
               <<"\"Q_criterion\" \\ \n"
	       <<"\"Vorticity\" \\ \n"
               <<"\"Fuel_rrate\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
      
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */

   
   
   for (int k = Grid.KNl ; k <= Grid.KNu ; ++k) {
      for (int j = Grid.JNl ; j <= Grid.JNu ; ++j) {
         for (int i = Grid.INl ; i <= Grid.INu ; ++i) {
	    W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X << W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T()
		     << " " << Q_criterion_n(*this, i, j, k)
		     << " " << vorticity_n(*this, i, j, k)
	             << " " << -W_node.Sw(W_node.React.reactset_flag).rhospec[0].c 
		     << "\n"; 
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
void Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
Output_Cells_Tecplot(Input_Parameters<LES3DTF_pState, 
                                      LES3DTF_cState> &IPs,
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
               << "\"k\" \\ \n";
      //n species mass fractions names
      for (int i =0; i<W[0][0][0].ns; ++i) {
         Out_File <<"\"c"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      } /* endif */
     
      Out_File <<"\"WF\" \\ \n" 
	       <<"\"TF\" \\ \n"
	       <<"\"Iso_c=0.5\" \\ \n";

      Out_File <<"\"T\" \\ \n"
               <<"\"R\" \\ \n";
      
      Out_File  << "\"mu_t\" \\ \n";

      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.ICu - Grid.ICl + 2*Grid.Nghost + 1 << " \\ \n"
               << "J = " << Grid.JCu - Grid.JCl + 2*Grid.Nghost + 1 << " \\ \n"
               << "K = " << Grid.KCu - Grid.KCl + 2*Grid.Nghost + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.ICu - Grid.ICl + 2*Grid.Nghost + 1 << " \\ \n"
               << "J = " << Grid.JCu - Grid.JCl + 2*Grid.Nghost + 1 << " \\ \n"
               << "K = " << Grid.KCu - Grid.KCl + 2*Grid.Nghost + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";
      
   } /* endif */

   for (int k = Grid.KCl-Grid.Nghost; k <= Grid.KCu+Grid.Nghost; ++k) {
      for (int j  = Grid.JCl-Grid.Nghost; j <= Grid.JCu+Grid.Nghost; ++j ) {
         for (int i = Grid.ICl-Grid.Nghost; i <= Grid.ICu+Grid.Nghost; ++i ) {
            
            Out_File << " "  <<  Grid.Cell[i][j][k].Xc
                     <<  W[i][j][k];
            Out_File.setf(ios::scientific);
            Out_File << " " <<W[i][j][k].T() 
                     << " " <<W[i][j][k].Rtot()<<"\n ";

            Out_File << " " <<W[i][j][k].mu_t(dWdx[i][j][k],
                                              dWdy[i][j][k],
                                              dWdz[i][j][k],
                                              Flow_Type,
                                              Grid.Cell[i][j][k].V) << "\n";
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
void Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
Output_Nodes_Tecplot(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs,
                     const int Number_of_Time_Steps,
                     const double &Time,
                     const int Block_Number,
                     const int Output_Title,
                     ostream &Out_File) {

   LES3DTF_pState W_node;

   /* Ensure boundary conditions are updated before
      evaluating solution at the nodes. */
   
   BCs(IPs);

   /* Output node solution data. */
   
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
               << "\"k\" \\ \n";
      //n species mass fractions names
      for(int i =0 ;i<W[0][0][0].ns ; ++i){
         Out_File <<"\"c_"<<W[0][0][0].specdata[i].Speciesname()<<"\" \\ \n";
      }
      Out_File <<"\"WF\" \\ \n" 
	       <<"\"TF\" \\ \n"
	       <<"\"Iso_c=0.5\" \\ \n";

      Out_File <<"\"T\" \\ \n";
      
      Out_File<< "ZONE T =  \"Block Number = " << Block_Number
              << "\" \\ \n"
              << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
              << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
              << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
              << "DATAPACKING = POINT \n";
   } else {
      Out_File << "ZONE T =  \"Block Number = " << Block_Number
               << "\" \\ \n"
               << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
               << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
               << "K = " << Grid.KNu - Grid.KNl + 1 << " \\ \n"
               << "DATAPACKING = POINT \n";              
   } /* endif */
   
   for (int k  = Grid.KNl; k <= Grid.KNu; ++k) {
      for (int j  = Grid.JNl; j <= Grid.JNu; ++j) {
         for (int i = Grid.INl; i <= Grid.INu; ++i) {
            W_node = Wn(i, j, k);
            Out_File << " "  << Grid.Node[i][j][k].X <<W_node;
            Out_File.setf(ios::scientific);
            Out_File << " " << W_node.T()<< "\n"; 

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
int Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
ICs(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs){
   
  LES3DTF_pState Wl, Wr;

  if (IPs.Wo.React.reactset_flag == CH4_1STEP) {
    //set to phi=1.0 values from CHEMKIN
    Wl.v.x = 0.3837;//0.4101;	
    Wr.v.x = 2.8582;//3.103; 

    //phi = 1.0
    Wl.spec[0] = 0.05518;       //CH4
    Wl.spec[1] = 0.22015;     //O2
    Wl.spec[2] = 0.0;     //CO2
    Wl.spec[3] = 0.0;     //H2O 
    Wl.spec[4] = 0.72467;
    Wl.p = 101325.0;
    Wl.rho = 1.13; //2234
    Wl.flame.TF = 1.0;
    Wl.flame.WF =1.0;

    Wr.spec[0] = ZERO;       //CH4
    Wr.spec[1] = 0.0000;     //O2
    Wr.spec[2] = 0.1511;     //CO2
    Wr.spec[3] = 0.1242;     //H2O 
    Wr.spec[4] = 0.72467;
    Wr.p = 101325.0;
    Wr.rho = Wr.p/(Wr.Rtot()*2230); //2234  2320
    Wr.flame.TF = 1.0;
    Wr.flame.WF =1.0;
 
  } else if (IPs.Wo.React.reactset_flag == CH4_2STEP){ 

    //phi = 1.0
    Wl.v.x = 0.4101;
    Wr.v.x = 3.103;
    Wr.spec[0] = ZERO;       //CH4
    Wr.spec[1] = 0.0055;      //O2
    Wr.spec[2] = 0.1368;     //CO2
    Wr.spec[3] = 0.1242;     //H2O 
    Wr.spec[4] = 0.0088;     //CO

    Wr.rho = Wr.p/(Wr.Rtot()*2250);
  }

  
  // begin switch
  switch(IPs.i_ICs) {

  case IC_UNIFORM :
    // Set the solution state to the initial state Wo[0].
    for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
      for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	  W[i][j][k] = IPs.Wo;
	  U[i][j][k] = W[i][j][k].U();
	} /* endfor */
      } /* endfor */
    } /* endfor */
    break;

  case IC_TURBULENT_BOX :
    for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
      for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	  W[i][j][k] = IPs.Wo;
	  W[i][j][k].rho = IPs.Wo.p/(IPs.Wo.Rtot()*298.15);
	  W[i][j][k].v.zero();
	  W[i][j][k].k = 0.0;
	  U[i][j][k] = W[i][j][k].U();
	} /* endfor */
      } /* endfor */
    } /* endfor */
    break;
      
  case IC_TURBULENT_PREMIXED_FLAME :     
     for (int k  = KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
        for (int j  = JCl- Nghost ; j <=  JCu+ Nghost ; ++j) {
           for (int i = ICl- Nghost ; i <=  ICu+ Nghost ; ++i) {
	     double xx = Grid.Cell[i][j][k].Xc.x;
	     double tau = Wr.T()/Wl.T() - ONE;
	     double C = 0.5*(ONE + erf(SQRT_PI*xx/(THREE*IPs.Wo._laminar_flame_thickness)));  // (erf(xx*2000.0)+1.0)/2.0;
	     if (C>0.02 && C<0.98) {
	       W[i][j][k].flame.TF = IPs.Wo._TFactor;  //5.0;
	     } else {
	       W[i][j][k].flame.TF = 1.0;
	     }
	     W[i][j][k].flame.WF =1.0;

	     double Yf_u = 0.05518;
	     double Yf_b = 0.0;
	     W[i][j][k].spec[0].c = Yf_b*C + Yf_u*(ONE-C);  
	     W[i][j][k].premixed_mfrac();

	     W[i][j][k].p = 101325.0;
	     W[i][j][k].rho = Wl.rho*Wl.Rtot()/W[i][j][k].Rtot()/(1.0+tau*C);
	     W[i][j][k].v.x = Wl.rho*0.3837/W[i][j][k].rho;
	     W[i][j][k].k = 0.00;
	     
	     U[i][j][k] = W[i][j][k].U();
	     
	   } /* endfor */
	} /* endfor */
     } /* endfor */
     break;

   case IC_TURBULENT_BUNSEN_BOX :
     {
     double C, xx, yy, zz, tau, Sl, slot_width;
     double Yf_u = 0.05518;
     double Yf_b = 0.0;
     for (int k  = KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
        for (int j  = JCl- Nghost ; j <=  JCu+ Nghost ; ++j) {
           for (int i = ICl- Nghost ; i <=  ICu+ Nghost ; ++i) {
	      xx = Grid.Cell[i][j][k].Xc.x;
	      yy = Grid.Cell[i][j][k].Xc.y;
	      zz = Grid.Cell[i][j][k].Xc.z;
	      tau = Wr.T()/Wl.T() - ONE;
	      Sl = IPs.Wo._laminar_flame_speed;
	      slot_width = 0.025;
	      
	      if (zz <= IPs.Fresh_Gas_Height) {
		if (yy < 0.0) {
		  C = 0.5*(ONE + erf(SQRT_PI*(-yy-HALF*slot_width)/(0.7*IPs.Wo._TFactor*IPs.Wo._laminar_flame_thickness)));
		} else {
		  C = 0.5*(ONE + erf(SQRT_PI*(yy-HALF*slot_width)/(0.7*IPs.Wo._TFactor*IPs.Wo._laminar_flame_thickness)));
		}
	      } else {
		C = ONE;
	      }

	      if (zz <= IPs.Fresh_Gas_Height && fabs(yy) <= HALF*slot_width) {
		// fresh gas inflow
		W[i][j][k].v.z = IPs.Mean_Velocity.z;  //15.58*(ONE-sqr(xx/0.005));
	      } else {
		W[i][j][k].v.z = IPs.Mean_Velocity.z + Sl*(Wl.rho/Wr.rho - ONE);
	      }

	      W[i][j][k].flame.TF = IPs.Wo._TFactor;
	      W[i][j][k].flame.WF =1.0;

	      W[i][j][k].spec[0].c = Yf_b*C + Yf_u*(ONE-C);
	      W[i][j][k].premixed_mfrac();
              W[i][j][k].p = 101325.0;
	      W[i][j][k].rho = Wl.rho*Wl.Rtot()/W[i][j][k].Rtot()/(1.0+tau*C);
	      W[i][j][k].k = 0.0;
	      W[i][j][k].v.x = 0.0;
	      W[i][j][k].v.y = 0.0;

	      if (zz <= IPs.Fresh_Gas_Height) {
		if (yy < 0.0) {
		  W[i][j][k].v.y = -Sl*(Wl.rho/W[i][j][k].rho - ONE);
		} else {
		  W[i][j][k].v.y = Sl*(Wl.rho/W[i][j][k].rho - ONE);
		}
	      }

	      if (fabs(yy) <= HALF*slot_width) {
		W[i][j][k].v.z = IPs.Mean_Velocity.z + Sl*(Wl.rho/W[i][j][k].rho - ONE); 
	      }

              U[i][j][k] = W[i][j][k].U();
	      
	   } /* endfor */
	} /* endfor */
     } /* endfor */
     }
     break;  

   case IC_TURBULENT_BUNSEN_FLAME :
     for (int k  = KCl- Nghost ; k <=  KCu+ Nghost ; ++k) {
        for (int j  = JCl- Nghost ; j <=  JCu+ Nghost ; ++j) {
           for (int i = ICl- Nghost ; i <=  ICu+ Nghost ; ++i) {
	     double rr = sqrt(sqr(Grid.Cell[i][j][k].Xc.x)+sqr(Grid.Cell[i][j][k].Xc.y));
	     double tau = Wr.T()/Wl.T() - ONE;
	     double C;

	     if (Grid.Cell[i][j][k].Xc.z<=0.035) {
	       C = 0.5*(ONE + erf(SQRT_PI*(rr-0.0056)/(THREE*0.44E-3)));
	       //C = ZERO;
	       if (rr<=0.0056){
		 W[i][j][k].v.z = 0.0; //15.58*(ONE-sqr(rr/0.0056));
	       }else{
		 W[i][j][k].v.z = 0.0;
		 W[i][j][k].v.x = 0.0;  //0.3837/sqrt(TWO);
	         W[i][j][k].v.y = 0.0;  //0.3837/sqrt(TWO);
	       } /* endif */
	     }else{
	       C = 1.0;
	       //C = ZERO;
	       W[i][j][k].v.z = 0.0;
	     } /* endif */

	     if (C>0.02 && C<0.98) {
	       W[i][j][k].flame.TF = 5.0;
	     } else {
	       W[i][j][k].flame.TF = 1.0;
	     }
	     //W[i][j][k].flame.TF = 5.0;
	     W[i][j][k].flame.WF =1.0;

	     double Yf_u = 0.05518;
	     double Yf_b = 0.0;
	     W[i][j][k].spec[0].c = Yf_b*C + Yf_u*(ONE-C);
	     W[i][j][k].premixed_mfrac();

	     W[i][j][k].p = 101325.0;
	     W[i][j][k].rho = 1.13*Wl.Rtot()/W[i][j][k].Rtot()/(1.0+tau*C);
	     W[i][j][k].v.x = 0.0; //1.13*0.3837/(W[i][j][k].rho*sqrt(TWO));
	     W[i][j][k].v.y = 0.0; //1.13*0.3837/(W[i][j][k].rho*sqrt(TWO));
	     W[i][j][k].k = 0.0;

	     U[i][j][k] = W[i][j][k].U();

	   } /* endfor */
	} /* endfor */
     } /* endfor */

     break;  
      
   } //end of switch
   
   /* Set default values for the boundary conditions
      reference states. */

   for (int k = KCl-Nghost ; k<= KCu+Nghost; ++k) {
      for (int j = JCl-Nghost ; j<= JCu+Nghost; ++j){
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
         } /* endif */
      } /* endfor */ 
   } /* endfor */
    
   for (int  k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
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
         } /* endif */
      } /* endfor */
   } /* endfor */

   for (int  j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
         if ((j >= JCl && j <= JCu) && (i >= ICl && i <= ICu)) {
            WoT[i][j] = W[i][j][KCu];
            WoB[i][j] = W[i][j][KCl];
         }  else if (i < ICl &&  j< JCl) {
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
         } /* endif */
      } /* endfor */
   } /* endfor */
      
   return (0);
    
}

/********************************************************
 * Routine: ICs_Specializations                         *
 *                                                      *
 * Apply initial conditions for the k-equation model    *
 *                                                      *
 ********************************************************/
template<>
int Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
ICs_Specializations(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs){

   /* Determine initial values for the subfilter scale
      turbulent kinetic energy. */
    
   Linear_Reconstruction_LeastSquares(IPs.i_Limiter);
   for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
      for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
          for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	    W[i][j][k].k = 0.005*sqr(W[i][j][k].filter_width()*
                               W[i][j][k].abs_strain_rate(dWdx[i][j][k],
                                                          dWdy[i][j][k],
                                                          dWdz[i][j][k]));
	    if (W[i][j][k].k < MICRO) W[i][j][k].k = ZERO; 
	    U[i][j][k] = W[i][j][k].U();
   	    } /* endfor */ 	  
   	} /* endfor */
   } /* endfor */

   /* Set default values for the boundary conditions
      reference states. */

   for (int k = KCl-Nghost ; k<= KCu+Nghost; ++k) {
      for (int j = JCl-Nghost ; j<= JCu+Nghost; ++j){
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
         } /* endif */
      } /* endfor */ 
   } /* endfor */
    
   for (int  k = KCl-Nghost ; k <= KCu+Nghost ; ++k) {
      for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
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
         } /* endif */
      } /* endfor */
   } /* endfor */

   for (int  j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
      for (int  i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
         if ((j >= JCl && j <= JCu) && (i >= ICl && i <= ICu)) {
            WoT[i][j] = W[i][j][KCu];
            WoB[i][j] = W[i][j][KCl];
         }  else if (i < ICl &&  j< JCl) {
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
         } /* endif */
      } /* endfor */
   } /* endfor */
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
void Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
BCs(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs) {

   double dpdx, dpdy, dpdz;
   Vector3D dX;
   Vector3D MOVING_WALL_VELOCITY = IPs.Moving_Wall_Velocity;

   for (int k = KCl- Nghost; k <=  KCu+ Nghost; ++k) {
      for (int j = JCl- Nghost; j <= JCu+ Nghost; ++j) {
         // Prescribe West boundary conditions.
         switch( Grid.BCtypeW[j][k]) {
           case BC_NONE :
             break;
            
           case BC_REFLECTION :
             W[ICl-1][j][k] = LES3DTF_pState::Reflect(W[ICl][j][k],Grid.nfaceW(ICl,j,k));
             U[ICl-1][j][k] = W[ICl-1][j][k].U();
             W[ICl-2][j][k] = LES3DTF_pState::Reflect(W[ICl+1][j][k],Grid.nfaceW(ICl,j,k));
             U[ICl-2][j][k] = W[ICl-2][j][k].U();
             break;

           case BC_NO_SLIP :
             W[ICl-1][j][k] = LES3DTF_pState::NoSlip(W[ICl][j][k],
                                                      WoW[j][k],
                                                      Grid.nfaceW(ICl,j,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[ICl-1][j][k] = W[ICl-1][j][k].U();
             W[ICl-2][j][k] = LES3DTF_pState::NoSlip(W[ICl+1][j][k],
                                                      WoW[j][k],
                                                      Grid.nfaceW(ICl,j,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[ICl-2][j][k] = W[ICl-2][j][k].U();
             break;

           case BC_MOVING_WALL :
             W[ICl-1][j][k] = LES3DTF_pState::MovingWall(W[ICl][j][k],
                                                          WoW[j][k],
                                                          Grid.nfaceW(ICl,j,k),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[ICl-1][j][k] = W[ICl-1][j][k].U();
             W[ICl-2][j][k] = LES3DTF_pState::MovingWall(W[ICl+1][j][k], 
                                                          WoW[j][k], 
                                                          Grid.nfaceW(ICl,j,k),
                                                          MOVING_WALL_VELOCITY, 
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[ICl-2][j][k] = W[ICl-2][j][k].U();
             break;

           case BC_FIXED_PRESSURE :
             W[ICl-1][j][k] = W[ICl][j][k];
             W[ICl-1][j][k].p = WoW[j][k].p;
             U[ICl-1][j][k] = W[ICl-1][j][k].U();
             W[ICl-2][j][k] = W[ICl][j][k] ;
             W[ICl-2][j][k].p = WoW[j][k].p;
             U[ICl-2][j][k] = W[ICl-2][j][k].U();
             break;

           case BC_CHANNEL_INFLOW:
             dpdx = IPs.Pressure_Gradient.x; 
             //for turbulent channel flow
             // p linearly varys based on constant pressure gradient 
             dX = Grid.Cell[ICl-1][j][k].Xc - Grid.Cell[ICl][j][k].Xc; 
             W[ICl-1][j][k] = WoW[j][k];
             W[ICl-1][j][k].v.x = W[ICl][j][k].v.x;
             W[ICl-1][j][k].p = WoW[j][k].p - dpdx*dX.x ;
             U[ICl-1][j][k] = W[ICl-1][j][k].U( );
             dX = Grid.Cell[ICl-2][j][k].Xc - Grid.Cell[ICl][j][k].Xc;
             W[ICl-2][j][k] = WoW[j][k];
             W[ICl-2][j][k].v.x = W[ICl][j][k].v.x;
             W[ICl-2][j][k].p = WoW[j][k].p - dpdx*dX.x;
             //W[ICl-2][j][k].p = (WoW[j][k].p - 2.0/3.0*W[ICl-2][j][k].rho*W[ICl-2][j][k].k) - dpdx*dX.x;
             U[ICl-2][j][k] = W[ICl-2][j][k].U( );
             break;
            
           case BC_PERIODIC :  
             W[ICl-1][j][k] = W[ICu-1][j][k];
             U[ICl-1][j][k] = U[ICu-1][j][k];
             W[ICl-2][j][k] = W[ICu-2][j][k];
             U[ICl-2][j][k] = U[ICu-2][j][k];
             break;

           case BC_INFLOW_SUBSONIC :
           // all fixed except pressure which is constant extrapolation
	     W[ICl-1][j][k] = WoW[j][k];
             W[ICl-1][j][k].p = W[ICl][j][k].p;
	     U[ICl-1][j][k] = W[ICl-1][j][k].U();
             W[ICl-2][j][k] = WoW[j][k];
             W[ICl-2][j][k].p = W[ICl][j][k].p;
	     U[ICl-2][j][k] = W[ICl-2][j][k].U();
	     break;

           case BC_OUTFLOW_SUBSONIC :
	   // all constant extrapolation except pressure which is fixed.
	     W[ICl-1][j][k] = W[ICl][j][k]; 
	     W[ICl-1][j][k].p = WoW[j][k].p;
	     U[ICl-1][j][k] = W[ICl-1][j][k].U();
	     
	     W[ICl-2][j][k].Copy( W[ICl-1][j][k] );  
	     U[ICl-2][j][k].Copy( U[ICl-1][j][k] );
	     break;

	   case BC_INFLOW_TURBULENCE :
	     W[ICl-1][j][k] = WoW[j][k];
	     W[ICl-1][j][k].v = W[ICl][j][k].v;
	     W[ICl-1][j][k].p = W[ICl][j][k].p;
	     W[ICl-1][j][k].k = W[ICl][j][k].k;
	     U[ICl-1][j][k] = W[ICl-1][j][k].U();

             W[ICl-2][j][k] = WoW[j][k];
	     W[ICl-2][j][k].v = W[ICl][j][k].v;
	     W[ICl-2][j][k].p = W[ICl][j][k].p;
	     W[ICl-2][j][k].k = W[ICl][j][k].k;
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
         switch( Grid.BCtypeE[j][k]) {
           case BC_NONE :
             break;

           case BC_REFLECTION :
             W[ICu+1][j][k] = LES3DTF_pState::Reflect(W[ICu][j][k],Grid.nfaceE(ICu,j,k));
             U[ICu+1][j][k] = W[ICu+1][j][k].U();
             W[ICu+2][j][k] = LES3DTF_pState::Reflect(W[ICu-1][j][k],Grid.nfaceE(ICu,j,k));
             U[ICu+2][j][k] = W[ICu+2][j][k].U();
             break;

           case BC_NO_SLIP :
             W[ICu+1][j][k] = LES3DTF_pState::NoSlip(W[ICu][j][k],
                                                      WoE[j][k],
                                                      Grid.nfaceE(ICu,j,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[ICu+1][j][k] = W[ICu+1][j][k].U();
             W[ICu+2][j][k] = LES3DTF_pState::NoSlip(W[ICu-1][j][k],
                                                      WoE[j][k],
                                                      Grid.nfaceE(ICu,j,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[ICu+2][j][k] = W[ICu+2][j][k].U();
             break;

           case BC_MOVING_WALL :
             W[ICu+1][j][k] = LES3DTF_pState::MovingWall(W[ICu][j][k],
                                                          WoE[j][k],
                                                          Grid.nfaceE(ICu,j,k),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[ICu+1][j][k] = W[ICu+1][j][k].U();
             W[ICu+2][j][k] = LES3DTF_pState::MovingWall(W[ICu-1][j][k], 
                                                          WoW[j][k], 
                                                          Grid.nfaceW(ICu,j,k),
                                                          MOVING_WALL_VELOCITY, 
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[ICu+2][j][k] = W[ICu+2][j][k].U();
             break;

           case BC_FIXED_PRESSURE :
             W[ICu+1][j][k] = W[ICu][j][k];
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
             U[ICu+2][j][k] = W[ICu+2][j][k].U( );
             break;
            
           case BC_PERIODIC :
             W[ICu+1][j][k] = W[ICl+1][j][k];
             U[ICu+1][j][k] = U[ICl+1][j][k];
             W[ICu+2][j][k] = W[ICl+2][j][k];
             U[ICu+2][j][k] = U[ICl+2][j][k];
             break;

 	   case BC_INFLOW_SUBSONIC :
             // all fixed except pressure which is constant extrapolation
	     W[ICu+1][j][k] = WoE[j][k];
             W[ICu+1][j][k].p = W[ICu][j][k].p;
	     U[ICu+1][j][k] = W[ICu+1][j][k].U();
 	     W[ICu+2][j][k] = WoE[j][k];
             W[ICu+2][j][k].p = W[ICu][j][k].p;
	     U[ICu+2][j][k] = W[ICu+2][j][k].U();
	     break;

	   case BC_OUTFLOW_SUBSONIC :
	     // all constant extrapolation except pressure which is fixed.
	     W[ICu+1][j][k] = W[ICu][j][k]; 
	     W[ICu+1][j][k].p = WoE[j][k].p;
	     U[ICu+1][j][k] = W[ICu+1][j][k].U();
	     
	     W[ICu+2][j][k].Copy( W[ICu+1][j][k] );
	     U[ICu+2][j][k].Copy( U[ICu+1][j][k] );
	     break;
            
	   case BC_INFLOW_TURBULENCE :
	     W[ICu+1][j][k] = WoE[j][k];
	     W[ICu+1][j][k].v = W[ICu][j][k].v;
	     W[ICu+1][j][k].p = W[ICu][j][k].p;
	     W[ICu+1][j][k].k = W[ICu][j][k].k;
	     U[ICu+1][j][k] = W[ICu+1][j][k].U();

 	     W[ICu+2][j][k] = WoE[j][k];
	     W[ICu+2][j][k].v = W[ICu][j][k].v;
	     W[ICu+2][j][k].p = W[ICu][j][k].p;
	     W[ICu+2][j][k].k = W[ICu][j][k].k;
	     U[ICu+2][j][k] = W[ICu+2][j][k].U();
	     break;

           case BC_CONSTANT_EXTRAPOLATION :
	   default :
             W[ICu+1][j][k] = W[ICu][j][k];
             U[ICu+1][j][k] = W[ICu+1][j][k].U();
             W[ICu+2][j][k] = W[ICu][j][k];
             U[ICu+2][j][k] = W[ICu+2][j][k].U();
             break;

         }/* endswitch */

      } /* endfor */
   } /* endfor */
       
   for (int k = KCl - Nghost; k <= KCu + Nghost; ++k) {
      for (int i = ICl - Nghost; i <= ICu + Nghost; ++i) {
         // Prescribe North boundary conditions.
         switch( Grid.BCtypeN[i][k]) {
           case BC_NONE :
             break;
            
           case BC_REFLECTION :
             W[i][JCu+1][k] = LES3DTF_pState::Reflect(W[i][JCu][k],Grid.nfaceN(i,JCu,k));
             U[i][JCu+1][k] = W[i][ JCu+1][k].U();
             W[i][JCu+2][k] = LES3DTF_pState::Reflect(W[i][JCu-1][k],Grid.nfaceN(i,JCu,k));
             U[i][JCu+2][k] = W[i][JCu+2][k].U();
             break;

           case BC_NO_SLIP :
             W[i][JCu+1][k] = LES3DTF_pState::NoSlip(W[i][JCu][k],
                                                      WoN[i][k],
                                                      Grid.nfaceN(i,JCu,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][JCu+1][k] = W[i][JCu+1][k].U();
             W[i][JCu+2][k] = LES3DTF_pState::NoSlip(W[i][JCu-1][k], 
                                                      WoN[i][k],
                                                      Grid.nfaceN(i,JCu,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][JCu+2][k] = W[i][JCu+2][k].U();
             break;
            
           case BC_MOVING_WALL :
             W[i][JCu+1][k] = LES3DTF_pState::MovingWall(W[i][JCu][k], 
                                                          WoN[i][k],
                                                          Grid.nfaceN(i,JCu,k),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][JCu+1][k] = W[i][JCu+1][k].U();
             W[i][JCu+2][k] = LES3DTF_pState::MovingWall(W[i][JCu-1][k], 
                                                          WoN[i][k],
                                                          Grid.nfaceN(i,JCu,k),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][JCu+2][k] = W[i][JCu+2][k].U();
             break;

           case BC_FIXED_PRESSURE :
             W[i][JCu+1][k] = W[i][JCu][k];
             W[i][JCu+1][k].p = WoN[i][k].p;
             U[i][JCu+1][k] = W[i][JCu+1][k].U();
             W[i][JCu+2][k] = W[i][JCu-1][k];
             W[i][JCu+2][k].p = WoN[i][k].p;
             U[i][JCu+2][k] = W[i][JCu+2][k].U();
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

           case BC_INFLOW_SUBSONIC :
	     // all fixed except p which is constant extrapolation
	     W[i][JCu+1][k] = WoN[i][k];
	     W[i][JCu+1][k].p = W[i][JCu][k].p;
	     U[i][JCu+1][k] = W[i][JCu+1][k].U();
 	     W[i][JCu+2][k] = WoN[i][k];
	     W[i][JCu+2][k].p = W[i][JCu][k].p;
	     U[i][JCu+2][k] = W[i][JCu+2][k].U();
	     break;

           case BC_OUTFLOW_SUBSONIC :
	     // all constant extrapolation except pressure which is fixed.
	     W[i][JCu+1][k] = W[i][JCu][k];
	     W[i][JCu+1][k].p = WoN[i][k].p;
	     U[i][JCu+1][k] = W[i][JCu+1][k].U();
	     
	     W[i][JCu+2][k].Copy( W[i][JCu+1][k] );
	     U[i][JCu+2][k].Copy( U[i][JCu+1][k] );
	     break;

	  case BC_INFLOW_TURBULENCE :
	     W[i][JCu+1][k] = WoN[i][k];
	     W[i][JCu+1][k].v = W[i][JCu][k].v;
	     W[i][JCu+1][k].p = W[i][JCu][k].p;
	     W[i][JCu+1][k].k = W[i][JCu][k].k;
	     U[i][JCu+1][k] = W[i][JCu+1][k].U();

 	     W[i][JCu+2][k] = WoN[i][k];
	     W[i][JCu+2][k].v = W[i][JCu][k].v;
	     W[i][JCu+2][k].p = W[i][JCu][k].p;
	     W[i][JCu+2][k].k = W[i][JCu][k].k;
	     U[i][JCu+2][k] = W[i][JCu+2][k].U();
	     break;

           case BC_CONSTANT_EXTRAPOLATION :
	   default :
             W[i][JCu+1][k] = W[i][JCu][k];
             U[i][JCu+1][k] = W[i][JCu+1][k].U();
             W[i][JCu+2][k] = W[i][JCu][k];
             U[i][JCu+2][k] = W[i][JCu+2][k].U();
             break;

         } /* endswitch */
    
         // Prescribe South boundary conditions.
         switch(Grid.BCtypeS[i][k]) {
           case BC_NONE :
             break;

           case BC_REFLECTION :
             W[i][JCl-1][k] = LES3DTF_pState::Reflect(W[i][JCl][k],Grid.nfaceS(i,JCl,k));
             U[i][JCl-1][k] = W[i][JCl-1][k].U();
             W[i][JCl-2][k] = LES3DTF_pState::Reflect(W[i][JCl+1][k],Grid.nfaceS(i,JCl,k));
             U[i][JCl-2][k] = W[i][JCl-2][k].U();
             break;

           case BC_NO_SLIP :
             W[i][JCl-1][k] = LES3DTF_pState::NoSlip(W[i][JCl][k], 
                                                      WoS[i][k], 
                                                      Grid.nfaceS(i,JCl,k),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][JCl-1][k] = W[i][JCl-1][k].U();
             W[i][JCl-2][k] = LES3DTF_pState::NoSlip(W[i][JCl+1][k], 
                                                      WoS[i][k],
                                                      Grid.nfaceS(i,JCl,k),   
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][JCl-2][k] = W[i][JCl-2][k].U();
             break;
 
           case BC_MOVING_WALL :
             W[i][JCl-1][k] = LES3DTF_pState::MovingWall(W[i][JCl][k], 
                                                          WoS[i][k], 
                                                          Grid.nfaceS(i,JCl,k),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient, 
                                                          FIXED_TEMPERATURE_WALL);
             U[i][JCl-1][k] = W[i][JCl-1][k].U();
             W[i][JCl-2][k] = LES3DTF_pState::MovingWall(W[i][JCl+1][k],
                                                          WoS[i][k],
                                                          Grid.nfaceS(i,JCl,k),
                                                          MOVING_WALL_VELOCITY, 
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][JCl-2][k] =  W[i][JCl-2][k].U();
             break;

           case BC_FIXED_PRESSURE :
             W[i][JCl-1][k] = W[i][JCl][k];
             W[i][JCl-1][k].p = WoS[i][k].p;
             U[i][JCl-1][k] =  W[i][JCl-1][k].U();
             W[i][JCl-2][k] = W[i][JCl+1][k];
             W[i][JCl-2][k].p = WoS[i][k].p;
             U[i][JCl-2][k] =  W[i][JCl-2][k].U();
	     break;

           case BC_CHANNEL_OUTFLOW:
             dpdy = IPs.Pressure_Gradient.y; 
             //for turbulent channel flow
             // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
             dX = Grid.Cell[i][JCl-1][k].Xc - Grid.Cell[i][JCl][k].Xc; 
             W[i][JCl-1][k] = WoS[i][k];
             W[i][JCl-1][k].p = WoS[i][k].p -dpdy*dX.y;
             //W[i][JCl-1][k].p = (WoS[i][k].p )-dpdy*dX.y;
             U[i][JCl-1][k] = W[i][JCl-1][k].U( );
             dX = Grid.Cell[i][JCl-2][k].Xc - Grid.Cell[i][JCl][k].Xc;
             W[i][JCl-2][k].p = WoS[i][k].p - dpdy*dX.y;
             //W[i][JCl-2][k].p = (WoS[i][k].p - 2.0/3.0*W[i][JCl-2][k].rho*W[i][JCl-2][k].k) - dpdy*dX.y;
             U[i][JCl-2][k] = W[i][JCl-2][k].U( );
             break;

           case BC_PERIODIC :
             W[i][JCl-1][k] = W[i][JCu-1][k];
             U[i][JCl-1][k] = U[i][JCu-1][k];
             W[i][JCl-2][k] = W[i][JCu-2][k];
             U[i][JCl-2][k] = U[i][JCu-2][k];
             break;

           case BC_INFLOW_SUBSONIC :
	     // all fixed except p which is constant extrapolation
	     W[i][JCl-1][k] = WoS[i][k];
	     W[i][JCl-1][k].p = W[i][JCl][k].p;
	     U[i][JCl-1][k] = W[i][JCl-1][k].U();
	     W[i][JCl-2][k] = WoS[i][k];
	     W[i][JCl-2][k].p = W[i][JCl][k].p;
	     U[i][JCl-2][k] = W[i][JCl-2][k].U();
	     break;

           case BC_OUTFLOW_SUBSONIC :
	     // all constant extrapolation except pressure which is fixed.
	     W[i][JCl-1][k] = W[i][JCl][k];
	     W[i][JCl-1][k].p = WoS[i][k].p;
	     U[i][JCl-1][k] = W[i][JCl-1][k].U();
	     
	     W[i][JCl-2][k].Copy( W[i][JCl-1][k] );
	     U[i][JCl-2][k].Copy( U[i][JCl-1][k] );
	     break;
            
	   case BC_INFLOW_TURBULENCE :
	     W[i][JCl-1][k] = WoS[i][k];
	     W[i][JCl-1][k].v = W[i][JCl][k].v;
	     W[i][JCl-1][k].p = W[i][JCl][k].p; 
	     W[i][JCl-1][k].k = W[i][JCl][k].k; 
	     U[i][JCl-1][k] = W[i][JCl-1][k].U();

	     W[i][JCl-2][k] = WoS[i][k];
	     W[i][JCl-2][k].v = W[i][JCl][k].v;
	     W[i][JCl-2][k].p = W[i][JCl][k].p; 
	     W[i][JCl-2][k].k = W[i][JCl][k].k; 
	     U[i][JCl-2][k] = W[i][JCl-2][k].U();
	     break;

           case BC_CONSTANT_EXTRAPOLATION :
	   default :
             W[i][JCl-1][k] = W[i][JCl][k];
             U[i][JCl-1][k] = W[i][JCl-1][k].U();
             W[i][JCl-2][k] = W[i][JCl][k];
             U[i][JCl-2][k] = W[i][JCl-2][k].U();
	     break;

         } /* endswitch */

      } /* endfor */
   } /* endfor */
 
   for (int j = JCl - Nghost; j <= JCu + Nghost; ++j) {
      for (int i = ICl - Nghost; i <= ICu + Nghost; ++i ) {
         // Prescribe Bottom boundary conditions.
         switch( Grid.BCtypeB[i][j]) {
           case BC_NONE :
             break;

           case BC_REFLECTION :
             W[i][j][KCl-1] = LES3DTF_pState::Reflect(W[i][j][KCl],Grid.nfaceBot(i,j,KCl));
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
             W[i][j][KCl-2] = LES3DTF_pState::Reflect(W[i][j][KCl+1],Grid.nfaceBot(i,j,KCl));
             U[i][j][KCl-2] = W[i][j][KCl-2].U();
             break;

           case BC_NO_SLIP :
             W[i][j][KCl-1] = LES3DTF_pState::NoSlip(W[i][j][KCl], 
                                                      WoB[i][j],
                                                      Grid.nfaceBot(i,j,KCl),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][j][KCl-1] = W[i][j][KCl-1].U(); 
             W[i][j][KCl-2] = LES3DTF_pState::NoSlip(W[i][j][KCl+1], 
                                                      WoB[i][j],
                                                      Grid.nfaceBot(i,j,KCl),
                                                      IPs.Pressure_Gradient,
                                                     FIXED_TEMPERATURE_WALL);
             U[i][j][KCl-2] = W[i][j][KCl-2].U();
             break;
            
           case BC_MOVING_WALL :
             W[i][j][KCl-1] = LES3DTF_pState::MovingWall(W[i][j][KCl], 
                                                          WoB[i][j],
                                                          Grid.nfaceBot(i,j,KCl),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
             W[i][j][KCl-2] = LES3DTF_pState::MovingWall(W[i][j][KCl+1], 
                                                          WoB[i][j],
                                                          Grid.nfaceBot(i,j ,KCl),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][j][KCl-2] =  W[i][j][KCl-2].U();
             break;

           case BC_FIXED_PRESSURE :
             W[i][j][KCl-1] = W[i][j][KCl];
             W[i][j][KCl-1].p = WoB[i][j].p;
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
             W[i][j][KCl-2] = W[i][j][KCl+1];
             W[i][j][KCl-2].p = WoB[i][j].p;
             U[i][j][KCl-2] = W[i][j][KCl-2].U();
             break;

           case BC_CHANNEL_INFLOW:
             dpdz = IPs.Pressure_Gradient.z; 
             //for turbulent channel flow
             // k and omega are constant extrapolation, p linearly varys based on constant pressure gradient 
             dX = Grid.Cell[i][j][KCl-1].Xc - Grid.Cell[i][j][KCl].Xc; 
             if((j==JCl-1 || j==JCl-2 || j== JCu+1 || j== JCu+2 ) && 
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
               U[i][j][KCl-1] = W[i][j][KCl-1].U( );
               dX = Grid.Cell[i][j][KCl-2].Xc - Grid.Cell[i][j][KCl].Xc;
               W[i][j][KCl-2] = WoB[i][j];
               W[i][j][KCl-2].v.z = W[i][j][KCl].v.z;
               W[i][j][KCl-2].p = WoB[i][j].p - dpdz*dX.z;
               // W[i][j][KCl-2].p = (WoB[i][j].p - 2.0/3.0*W[i][j][KCl-2].rho* W[i][j][KCl-2].k ) - dpdz*dX.z;
               U[i][j][KCl-2] = W[i][j][KCl-2].U( );
             } /* endif */
             break;

           case BC_PERIODIC :
             W[i][j][KCl-1] = W[i][j][KCu-1];
             U[i][j][KCl-1] = U[i][j][KCu-1];
             W[i][j][KCl-2] = W[i][j][KCu-2];
             U[i][j][KCl-2] = U[i][j][KCu-2];
             break;

           case BC_INFLOW_SUBSONIC :
	     // all fixed except p which is constant extrapolation
             W[i][j][KCl-1] = WoB[i][j]; 
             W[i][j][KCl-1].p = W[i][j][KCl].p;
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
	     
	     W[i][j][KCl-2].Copy( W[i][j][KCl-1] );
             U[i][j][KCl-2].Copy( U[i][j][KCl-1] );
             break;

           case BC_OUTFLOW_SUBSONIC :
	     // all constant extrapolation except pressure which is fixed.
             W[i][j][KCl-1] = W[i][j][KCl];
             W[i][j][KCl-1].p = WoB[i][j].p;
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
             W[i][j][KCl-2] = W[i][j][KCl];
             W[i][j][KCl-2].p = WoB[i][j].p;
             U[i][j][KCl-2] = W[i][j][KCl-2].U();
             break;

           case BC_INFLOW_TURBULENCE :
             W[i][j][KCl-1] = WoB[i][j];
	     W[i][j][KCl-1].v = W[i][j][KCl].v;
	     W[i][j][KCl-1].p = W[i][j][KCl].p; 
	     W[i][j][KCl-1].k = W[i][j][KCl].k; 
             U[i][j][KCl-1] = W[i][j][KCl-1].U();

	     W[i][j][KCl-2].Copy( W[i][j][KCl-1] ); 
             U[i][j][KCl-2].Copy( U[i][j][KCl-1] );
             break;

           case BC_CONSTANT_EXTRAPOLATION :
	   default :
             W[i][j][KCl-1] = W[i][j][KCl];
             U[i][j][KCl-1] = W[i][j][KCl-1].U();
             W[i][j][KCl-2] = W[i][j][KCl];
             U[i][j][KCl-2] = W[i][j][KCl-2].U();
             break;

         } /* endswitch */
       
         // Prescribe Top boundary conditions.
         switch( Grid.BCtypeT[i][j]) {
           case BC_NONE :
             break;
            
           case BC_REFLECTION :
             W[i][j][KCu+1] = LES3DTF_pState::Reflect(W[i][j][KCu],
                                                       Grid.nfaceTop(i,j,KCu));
             U[i][j][KCu+1] = W[i][j][KCu+1].U();
             W[i][j][KCu+2] = LES3DTF_pState::Reflect(W[i][j][KCu-1],
                                                       Grid.nfaceTop(i,j,KCu));
             U[i][j][KCu+2] = W[i][j][KCu+2].U();
             break;

           case BC_NO_SLIP :
             W[i][j][KCu+1] = LES3DTF_pState::NoSlip(W[i][j][KCu], 
                                                      WoT[i][j],
                                                      Grid.nfaceTop(i,j,KCu),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][j][KCu+1] = W[i][j][KCu+1].U();
             W[i][j][KCu+2] = LES3DTF_pState::NoSlip(W[i][j][KCu-1], 
                                                      WoT[i][j],
                                                      Grid.nfaceTop(i,j,KCu),
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
             U[i][j][KCu+2] = W[i][j][KCu+2].U();
             break;
            
           case BC_MOVING_WALL :
             W[i][j][KCu+1] = LES3DTF_pState::MovingWall(W[i][j][KCu], 
                                                          WoT[i][j],
                                                          Grid.nfaceTop(i,j,KCu),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][j][KCu+1] = W[i][j][KCu+1].U();
             W[i][j][KCu+2] = LES3DTF_pState::MovingWall(W[i][j][KCu-1], 
                                                          WoT[i][j],
                                                          Grid.nfaceTop(i,j,KCu),
                                                          MOVING_WALL_VELOCITY,
                                                          IPs.Pressure_Gradient,
                                                          FIXED_TEMPERATURE_WALL);
             U[i][j][KCu+2] = W[i][j][KCu+2].U();
             break;

           case BC_FIXED_PRESSURE :
             W[i][j][KCu+1] = W[i][j][KCu];
             W[i][j][KCu+1].p = WoT[i][j].p;
             U[i][j][KCu+1] =  W[i][j][KCu+1].U();
             W[i][j][KCu+2] =  W[i][j][KCu-1];
             W[i][j][KCu+2].p = WoT[i][j].p;
             U[i][j][KCu+2] =  W[i][j][KCu+2].U();
             break;

           case BC_CHANNEL_OUTFLOW:
             dpdz = IPs.Pressure_Gradient.z; 
             // all constant extrapolation except pressure specified which linearly varys if there is pressure gradient
             dX = Grid.Cell[i][j][KCu+1].Xc - Grid.Cell[i][j][KCu].Xc; 
             W[i][j][KCu+1] = W[i][j][KCu]; 
             W[i][j][KCu+1].p = (WoT[i][j].p)-dpdz*dX.z;
             U[i][j][KCu+1] = W[i][j][KCu+1].U( );
             dX = Grid.Cell[i][j][KCu+2].Xc - Grid.Cell[i][j][KCu].Xc; 
             W[i][j][KCu+2] = W[i][j][KCu];
             W[i][j][KCu+2].p = (WoT[i][j].p)-dpdz*dX.z; 	
             U[i][j][KCu+2] = W[i][j][KCu+2].U( );
             break;

           case BC_PERIODIC :
             W[i][j][KCu+1] = W[i][j][KCl+1];
             U[i][j][KCu+1] = U[i][j][KCl+1];
             W[i][j][KCu+2] = W[i][j][KCl+2];
             U[i][j][KCu+2] = U[i][j][KCl+2];
             break;

           case BC_INFLOW_SUBSONIC :
	     // all fixed except p which is constant extrapolation
             W[i][j][KCu+1] = WoT[i][j];
             W[i][j][KCu+1].p = W[i][j][KCu].p;
             U[i][j][KCu+1] = W[i][j][KCu+1].U();
             W[i][j][KCu+2] = WoT[i][j];
             W[i][j][KCu+2].p = W[i][j][KCu].p;
             U[i][j][KCu+2] = W[i][j][KCu+2].U();
             break;

           case BC_OUTFLOW_SUBSONIC :
	     // all constant extrapolation except pressure which is fixed.
             W[i][j][KCu+1] = W[i][j][KCu];
             W[i][j][KCu+1].p = WoT[i][j].p;
             U[i][j][KCu+1] = W[i][j][KCu+1].U();
	     
	     W[i][j][KCu+2].Copy( W[i][j][KCu+1] );
             U[i][j][KCu+2].Copy( U[i][j][KCu+1] );
             break;

	   case BC_INFLOW_TURBULENCE :
	     W[i][j][KCu+1] = WoT[i][j];
	     W[i][j][KCu+1].v = W[i][j][KCu].v; 
	     W[i][j][KCu+1].p = W[i][j][KCu].p; 
	     W[i][j][KCu+1].k = W[i][j][KCu].k; 
             U[i][j][KCu+1] = W[i][j][KCu+1].U();

             W[i][j][KCu+2] = WoT[i][j];
	     W[i][j][KCu+2].v = W[i][j][KCu].v; 
	     W[i][j][KCu+2].p = W[i][j][KCu].p; 
	     W[i][j][KCu+2].k = W[i][j][KCu].k; 
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
double Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
CFL(Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs) {
   
   double dtMin, d_i, d_j, d_k, v_i, v_j, v_k, a, dt_vis, nu, nu_t;
   double mr, aa_i, aa_j, aa_k;
   double length_n, delta_n, dTime;
   Vector3D V;
   
   dtMin = MILLION;
   
   for (int k = KCl- Nghost; k <= KCu + Nghost; ++k) {
      for (int j =  JCl- Nghost; j <= JCu + Nghost; ++j) {
         for (int i =  ICl- Nghost; i <= ICu+ Nghost; ++i) {
            if (i <  ICl || i >  ICu ||
                j <  JCl || j >  JCu || 
                k <  KCl || k >  KCu) {
               dt[i][j][k] = ZERO;
            } else {
	       V = W[i][j][k].rhov()/W[i][j][k].rho;

               d_i = TWO*(Grid.Cell[i][j][k].V/
                          (Grid.AfaceE(i, j, k)+ Grid.AfaceW(i, j, k)));
               d_j = TWO*( Grid.Cell[i][j][k].V/
                           (Grid.AfaceN(i, j, k)+ Grid.AfaceS(i, j, k)));
               d_k = TWO*( Grid.Cell[i][j][k].V/
                           (Grid.AfaceTop(i, j, k)+ Grid.AfaceBot(i, j, k)));
               v_i = HALF*(V * (Grid.nfaceE(i, j, k)- Grid.nfaceW(i, j, k)));
               v_j = HALF*(V * (Grid.nfaceN(i, j, k)- Grid.nfaceS(i, j, k)));
               v_k = HALF*(V * (Grid.nfaceTop(i, j, k)- Grid.nfaceBot(i, j, k)));

  	       //no preconditioning
               if(IPs.Preconditioning == 0){
                 a =  W[i][j][k].a();
                 dt[i][j][k] = min(min(d_i/(a+fabs(v_i)), d_j/(a+fabs(v_j))),
                                   d_k/(a+fabs(v_k)));
               
	       //Low Mach Number Preconditioning
 	       } else if (IPs.Preconditioning == 1) {
		 length_n = max(max(IPs.Grid_IP.Box_Length, 
                                  IPs.Grid_IP.Box_Width), 
                                  IPs.Grid_IP.Box_Height);
		 delta_n = min(min(fabs(d_i),fabs(d_j)),fabs(d_k));

	         dt[i][j][k] = min(min(d_i/W[i][j][k].u_plus_aprecon(fabs(v_i),
     		 					             delta_n,
                                                                     length_n,
                                                                     dTime),
	       			       d_j/W[i][j][k].u_plus_aprecon(fabs(v_j),
					    		             delta_n,
                                                                     length_n,
								     dTime)),
	       			       d_k/W[i][j][k].u_plus_aprecon(fabs(v_k),
					  		             delta_n,
                                                                     length_n,
                                                                     dTime));
	       } /* endif */
                     
               if (Flow_Type != FLOWTYPE_INVISCID) {  
                  nu = W[i][j][k].mu()/W[i][j][k].rho;

                  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K || 
		      Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {  
		     nu_t = W[i][j][k].mu_t(dWdx[i][j][k],
                                            dWdy[i][j][k],
                                            dWdz[i][j][k],
                                            Flow_Type,
					    Grid.Cell[i][j][k].V) / W[i][j][k].rho; 
                     nu = max(nu, nu_t);
		  }
                  
                  dt_vis = min(min((d_i*d_i)/nu, (d_j*d_j)/nu), (d_k*d_k)/nu)/3.0; 
                  dt[i][j][k]  = min(dt_vis, dt[i][j][k]);
		  
             /*
             int nn = NumVar();
             DenseMatrix dSwdU(nn,nn);
             dSwdU.zero();
             double max_diagonal = ONE;

             W[i][j][k].SemiImplicitSourceJacobi(dWdx[i][j][k],dWdy[i][j][k],dWdz[i][j][k],
                                                 d_dWdx_dW[i][j][k][0],d_dWdy_dW[i][j][k][0],
                                                 d_dWdz_dW[i][j][k][0],dSwdU,Grid.volume(i,j,k)); 

             if(IPs.Preconditioning == 1){
 	      DenseMatrix Pinv(nn,nn);
 	      Pinv.zero();
 	      W[i][j][k].Low_Mach_Number_Preconditioner_Inverse(Pinv, delta_n, length_n, dTime);	 
 	      dSwdU = Pinv*dSwdU;
             }
             for(int ii=0; ii<nn; ii++){
              max_diagonal = max(max_diagonal,fabs(dSwdU(ii,ii)));
              }
 	       dt[i][j][k] = min(HALF/max_diagonal, dt[i][j][k]);
		  */

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
int Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
dUdt_Multistage_Explicit(const int i_stage,
                         Input_Parameters<LES3DTF_pState, LES3DTF_cState> &IPs) {

   int i, j, k,  k_residual;
   double omega;
   Vector3D dX;
   
   LES3DTF_pState Wl, Wr;
   LES3DTF_cState Flux, Temp;

   LES3DTF_cState U_VACUUM;
   U_VACUUM.Vacuum();
   LES3DTF_pState W_VACUUM;
   W_VACUUM.Vacuum();

   double delta_n;

   Wl = W_VACUUM;
   Wr = W_VACUUM;
   Flux = U_VACUUM;
   Temp = U_VACUUM;

   // Evaluate the time step fraction and residual storage location for the stage. 
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

   /* Evaluate the time rate of change of the solution
      (i.e., the solution residuals) using a second-order
      limited upwind scheme with a variety of flux functions. */
   
   // Add i-direction (zeta-direction) fluxes.
   for ( k = KCl-1; k <= KCu+1; ++k ) {
      for ( j = JCl-1; j <= JCu+1; ++j ) {
         if ( i_stage == 1 ) {
            Uo[ICl-1][j][k] = U[ICl-1][j][k];
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
         } else {
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
         } /* endif */
         
         for ( i = ICl-1; i <= ICu; ++i ) {
            if ( i_stage == 1 ) {
               Uo[i+1][j][k] = U[i+1][j][k];
               dUdt[i+1][j][k][k_residual] = U_VACUUM;
            } else if (( j > JCl-1 && j < JCu+1 ) 
                       && (k > KCl-1 && k < KCu+1)) {
               switch(IPs.i_Time_Integration) {
                 case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:
                   //  dUdt[i+1][j][k][k_residual] = 
                   //  dUdt[i+1][j][k][k_residual];
                   break;

                 case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA:
		   if (IPs.N_Stage == 2) {
                     // dUdt[i+1][j][k_residual] = 
                     // dUdt[i+1][j][k_residual];
                   } else if (IPs.N_Stage == 4 && i_stage == 4) {
                      dUdt[i+1][j][k][k_residual] = dUdt[i+1][j][k][0] + 
                                                    TWO*dUdt[i+1][j][k][1] +
                                                    TWO*dUdt[i+1][j][k][2];
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

            if (( j > JCl-1 && j < JCu+1 ) && 
                ( k > KCl-1 && k < KCu+1 )) {

               /* Evaluate the cell interface i-direction fluxes. */
               if (i ==  ICl-1 && (Grid.BCtypeW[j][k] == BC_REFLECTION ||
                                   Grid.BCtypeW[j][k] == BC_NO_SLIP||
                                   Grid.BCtypeW[j][k] == BC_MOVING_WALL) ) {
                  dX = Grid.xfaceW(i+1,j,k)- Grid.Cell[i+1][j][k].Xc;
                  Wr = W[i+1][j][k] + 
                     ( phi[i+1][j][k]^dWdx[i+1][j][k])*dX.x +
                     ( phi[i+1][j][k]^dWdy[i+1][j][k])*dX.y +
                     ( phi[i+1][j][k]^dWdz[i+1][j][k])*dX.z;
                                 
                  if ( Grid.BCtypeW[j][k] == BC_REFLECTION) {
                     Wl = LES3DTF_pState::Reflect(Wr,Grid.nfaceW(i+1,j,k));
                  } /* endif */

                  if ( Grid.BCtypeW[j][k] == BC_NO_SLIP) {
                     Wl = LES3DTF_pState::NoSlip(Wr,WoW[j][k],
                                                  Grid.nfaceW(i+1,j,k), 
                                                  IPs.Pressure_Gradient,
                                                  FIXED_TEMPERATURE_WALL);
                                                         
                  } /* endif */

                  if ( Grid.BCtypeW[j][k] == BC_MOVING_WALL) {
                     Wl = LES3DTF_pState::MovingWall(Wr,WoW[j][k],
                                                      Grid.nfaceW(i+1,j,k),
                                                      IPs.Moving_Wall_Velocity,
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
                  } /* endif */

               } else if (i ==  ICu && 
                         (Grid.BCtypeE[j][k] == BC_REFLECTION||
                          Grid.BCtypeE[j][k] == BC_NO_SLIP||
                          Grid.BCtypeE[j][k] == BC_MOVING_WALL)) {
                  dX = Grid.xfaceE(i,j,k)- Grid.Cell[i][j][k].Xc;
                  Wl = W[i][j][k] + 
                     ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                     ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                     ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;

                  if ( Grid.BCtypeE[j][k] == BC_REFLECTION) {
                     Wr = LES3DTF_pState::Reflect(Wl,Grid.nfaceE(i,j,k));
                  } /* endif */

                  if ( Grid.BCtypeE[j][k] == BC_NO_SLIP) {
                     Wr = LES3DTF_pState::NoSlip(Wl,WoE[j][k],
                                                  Grid.nfaceE(i,j,k), 
                                                  IPs.Pressure_Gradient, 
                                                  FIXED_TEMPERATURE_WALL);
                  } /* endif */

                  if ( Grid.BCtypeE[j][k] == BC_MOVING_WALL) {
                     Wr = LES3DTF_pState::MovingWall(Wl,WoE[j][k],
                                                      Grid.nfaceE(i,j,k),  
                                                      IPs.Moving_Wall_Velocity,
                                                      IPs.Pressure_Gradient,
                                                      FIXED_TEMPERATURE_WALL);
                  } /* endif */

               } else { 
                  dX = Grid.xfaceE(i,j,k)- Grid.Cell[i][j][k].Xc;
                  Wl = W[i][j][k] + 
                       ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                       ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                       ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
                  dX = Grid.xfaceW(i+1,j,k)- Grid.Cell[i+1][j][k].Xc;
                  Wr = W[i+1][j][k] + 
                       ( phi[i+1][j][k]^ dWdx[i+1][j][k])*dX.x +
                       ( phi[i+1][j][k]^ dWdy[i+1][j][k])*dX.y +
                       ( phi[i+1][j][k]^ dWdz[i+1][j][k])*dX.z;

               } /* endif */
               
               // Spacing for Preconditioner 
	       if(IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
 	          delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
                                    TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                    TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	       } /* endif */

               // Evaluate inviscid numerical flux
               switch(IPs.i_Flux_Function) {
                case FLUX_FUNCTION_HLLE :
                  Flux = LES3DTF_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceE(i,j,k));
                  break;
                case FLUX_FUNCTION_ROE :
                  Flux = LES3DTF_pState::FluxRoe_n(Wl, Wr, Grid.nfaceE(i,j,k));
                  break;
   	        case FLUX_FUNCTION_AUSM_PLUS_UP :
 		  Flux = LES3DTF_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceE(i,j,k));
 	          break;
               } /* endswitch */

               // Add viscous flux in i direction
               Flux -=  LES3DTF_pState::FluxViscous_n(Wl, Wr,
                                                       W[i][j][k], W[i+1][j][k], 
                                                       dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                       dWdx[i+1][j][k], dWdy[i+1][j][k], dWdz[i+1][j][k],
                                                       Grid.nfaceE(i,j,k), Grid.Voe(i,j,k),  
                                                       Grid.delta_oe(i,j,k), Grid.Cell[i][j][k].V, 
	                                               Grid.Cell[i+1][j][k].V, Flow_Type);

               /* Evaluate cell-averaged solution changes. */
               dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                                            Flux*Grid.AfaceE(i,j,k)/Grid.Cell[i][j][k].V;
               dUdt[i+1][j][k][k_residual] += (IPs.CFL_Number*dt[i+1][j][k])*
                                              Flux*Grid.AfaceW(i+1,j,k)/Grid.Cell[i+1][j][k].V;

               /* Include source terms associated with modelled turbulence-chemistry interactions. */
	       dUdt[i][j][k][k_residual] += (IPs.CFL_Number*dt[i][j][k])*W[i][j][k].Sw(W[i][j][k].React.reactset_flag);
	       

               dUdt[i][j][k][k_residual] += (IPs.CFL_Number*dt[i][j][k])*
                                            LES3DTF_pState::Sturbchem(W[i][j][k],
                                                                      dWdx[i][j][k],
                                                                      dWdy[i][j][k], 
                                                                      dWdz[i][j][k],
                                                                      Flow_Type,
                                                                      Grid.Cell[i][j][k].V);

	       /* Include source terms associated with gravity */
	       //dUdt[i][j][k][k_residual] += (IPs.CFL_Number*dt[i][j][k])*W[i][j][k].Sg();


               /* Include physical time derivative for dual time stepping. */
//                if (IPs.Dual_Time_Stepping) {
//                   dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
//                                                W[i][j][k].S_dual_time_stepping(U[i][j][k], Ut[i][j][k], Uold[i][j][k], 
//                                                                                dTime, IPs.first_step);
//                } /* endif */

               /* Save west and east face boundary flux. */
               //    if (i ==  ICl-1) {
//                     FluxW[j] = -Flux* Grid.lfaceW(i+1, j);
//                 } else if (i ==  ICu) {
//                     FluxE[j] = Flux* Grid.lfaceE(i, j);
//                 } /* endif */ 
               

            } /* endif */
         } /* endfor */
          
         if (( j > JCl-1 && j < JCu+1 ) && 
             ( k > KCl-1 && k < KCu+1 )) {
            dUdt[ICl-1][j][k][k_residual] = U_VACUUM;
            dUdt[ICu+1][j][k][k_residual] = U_VACUUM;
         } /* endif */
      } /* endfor */
   } /* endfor */
 
   // Add j-direction (eta-direction) fluxes.
   for ( k = KCl; k <= KCu; ++k ) {
      for ( i = ICl; i <= ICu; ++i ) {
         for ( j = JCl-1; j <= JCu; ++j ) {

            /* Evaluate the cell interface j-direction fluxes. */
            if (j == JCl-1 && 
                ( Grid.BCtypeS[i][k] == BC_REFLECTION ||
                  Grid.BCtypeS[i][k] == BC_NO_SLIP||
                  Grid.BCtypeS[i][k] == BC_MOVING_WALL)) {
               dX = Grid.xfaceS(i,j+1,k)- Grid.Cell[i][j+1][k].Xc;
               Wr = W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y+
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;

               if ( Grid.BCtypeS[i][k] == BC_REFLECTION) {
                  Wl = LES3DTF_pState::Reflect(Wr,Grid.nfaceS(i,j+1,k));
	       } /* endif */

               if ( Grid.BCtypeS[i][k] == BC_NO_SLIP) {
                  Wl = LES3DTF_pState::NoSlip(Wr,WoS[i][k],
                                               Grid.nfaceS(i,j+1,k),
                                               IPs.Pressure_Gradient,
                                               FIXED_TEMPERATURE_WALL);
	       } /* endif */

               if ( Grid.BCtypeS[i][k] == BC_MOVING_WALL) {
                  Wl = LES3DTF_pState::MovingWall(Wr,WoS[i][k],
                                                   Grid.nfaceS(i,j+1,k),
                                                   IPs.Moving_Wall_Velocity,
                                                   IPs.Pressure_Gradient,
                                                   FIXED_TEMPERATURE_WALL);
	       } /* endif */

            } else if ( j ==  JCu && 
                      ( Grid.BCtypeN[i][k] == BC_REFLECTION ||
                        Grid.BCtypeN[i][k] == BC_NO_SLIP||
                        Grid.BCtypeN[i][k] == BC_MOVING_WALL)) {
               dX = Grid.xfaceN(i,j,k)- Grid.Cell[i][j][k].Xc;
               Wl = W[i][j][k] + 
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;

               if ( Grid.BCtypeN[i][k] == BC_REFLECTION) {
                  Wr = LES3DTF_pState::Reflect(Wl,Grid.nfaceN(i,j,k));
               } /* endif */

               if ( Grid.BCtypeN[i][k] == BC_NO_SLIP) {
                  Wr = LES3DTF_pState::NoSlip(Wl,WoN[i][k],
                                               Grid.nfaceN(i,j,k),
                                               IPs.Pressure_Gradient,
                                               FIXED_TEMPERATURE_WALL);
               } /* endif */

               if ( Grid.BCtypeN[i][k] == BC_MOVING_WALL) {
                  Wr = LES3DTF_pState::MovingWall(Wl,WoN[i][k],
                                                   Grid.nfaceN(i,j,k),
                                                   IPs.Moving_Wall_Velocity,
                                                   IPs.Pressure_Gradient,
                                                   FIXED_TEMPERATURE_WALL);
               } /* endif */

            } else {
               dX = Grid.xfaceN(i,j,k)- Grid.Cell[i][j][k].Xc;
               Wl = W[i][j][k] + 
                  ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                  ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                  ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               dX =  Grid.xfaceS(i,j+1,k)- Grid.Cell[i][j+1][k].Xc;
               Wr =  W[i][j+1][k] +
                  ( phi[i][j+1][k]^ dWdx[i][j+1][k])*dX.x +
                  ( phi[i][j+1][k]^ dWdy[i][j+1][k])*dX.y +
                  ( phi[i][j+1][k]^ dWdz[i][j+1][k])*dX.z;

	    } /* endif */

	    // Spacing for Preconditioner 
	    if (IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning){
 	       delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
                                 TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                 TWO*( Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	    } /* endif */

            // Evaluate inviscid numerical flux
            switch(IPs.i_Flux_Function) {
              case FLUX_FUNCTION_HLLE :
                Flux =  LES3DTF_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceN(i,j,k));
                break;
              case FLUX_FUNCTION_ROE :
                Flux =  LES3DTF_pState::FluxRoe_n(Wl, Wr, Grid.nfaceN(i,j,k));
                break;
              case FLUX_FUNCTION_AUSM_PLUS_UP :
                Flux =  LES3DTF_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceN(i,j,k));
	        break;
            } /* endswitch */

            // Add viscous flux in j direction
            Flux -= LES3DTF_pState::FluxViscous_n(Wl, Wr, 
                                                   W[i][j][k], W[i][j+1][k],
                                                   dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                   dWdx[i][j+1][k], dWdy[i][j+1][k], dWdz[i][j+1][k],
                                                   Grid.nfaceN(i,j,k), Grid.Von(i,j,k),
                                                   Grid.delta_on(i,j,k), Grid.Cell[i][j][k].V, 
                                                   Grid.Cell[i][j+1][k].V, Flow_Type);
            
            /* Evaluate cell-averaged solution changes. */
            dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                                         Flux*Grid.AfaceN(i,j,k)/Grid.Cell[i][j][k].V;
            dUdt[i][j+1][k][k_residual] += (IPs.CFL_Number*dt[i][j+1][k])*
                                           Flux*Grid.AfaceS(i,j+1,k)/Grid.Cell[i][j+1][k].V;

            /* Save south and north face boundary flux. */

//           if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */

	 } /* endfor */

         dUdt[i][JCl-1][k][k_residual] = U_VACUUM;
         dUdt[i][JCu+1][k][k_residual] = U_VACUUM;
      } /* endfor */
   } /* endfor */

   // Add k-direction (gamma-direction) fluxes.
   for ( i = ICl; i <= ICu; ++i ) {
      for ( j = JCl; j <= JCu; ++j ){
         for ( k = KCl-1; k <= KCu; ++k )  {

            /* Evaluate the cell interface j-direction fluxes. */
            if ( k == KCl-1 && 
               ( Grid.BCtypeB[i][j] == BC_REFLECTION  || 
                 Grid.BCtypeB[i][j] == BC_NO_SLIP ||
                 Grid.BCtypeB[i][j] == BC_MOVING_WALL)) {
                
               dX = Grid.xfaceBot(i,j,k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr = W[i][j][k+1] +
                  ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                  ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y+
                  ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;

               if ( Grid.BCtypeB[i][j] == BC_REFLECTION) {
                  Wl = LES3DTF_pState::Reflect(Wr,Grid.nfaceBot(i,j,k+1));
               } /* endif */

               if ( Grid.BCtypeB[i][j] == BC_NO_SLIP) {
                  Wl = LES3DTF_pState::NoSlip(Wr,WoB[i][j],
                                               Grid.nfaceBot(i,j,k+1), 
                                               IPs.Pressure_Gradient,
                                               FIXED_TEMPERATURE_WALL);
               } /* endif */

               if ( Grid.BCtypeB[i][j] == BC_MOVING_WALL) {
                  Wl = LES3DTF_pState::MovingWall(Wr,WoB[i][j],
                                                   Grid.nfaceBot(i,j,k+1),
                                                   IPs.Moving_Wall_Velocity,
                                                   IPs.Pressure_Gradient,
                                                   FIXED_TEMPERATURE_WALL);
               } /* endif */

            } else if ( k == KCu && 
                      ( Grid.BCtypeT[i][j] == BC_REFLECTION ||
                        Grid.BCtypeT[i][j] == BC_NO_SLIP ||
                        Grid.BCtypeT[i][j] == BC_MOVING_WALL)) {
                
                dX = Grid.xfaceTop(i,j,k)- Grid.Cell[i][j][k].Xc;
                Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x+
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y+
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;

                if ( Grid.BCtypeT[i][j] == BC_REFLECTION) {
                   Wr = LES3DTF_pState::Reflect(Wl,Grid.nfaceTop(i,j,k));
                } /* endif */

                if ( Grid.BCtypeT[i][j] == BC_NO_SLIP) {
                   Wr = LES3DTF_pState::NoSlip(Wl,WoT[i][j],
                                                Grid.nfaceTop(i,j,k), 
                                                IPs.Pressure_Gradient,
                                                FIXED_TEMPERATURE_WALL );
                } /* endif */

                if ( Grid.BCtypeT[i][j] == BC_MOVING_WALL) {
                   Wr = LES3DTF_pState::MovingWall(Wl,WoT[i][j],
                                                    Grid.nfaceTop(i,j,k),
                                                    IPs.Moving_Wall_Velocity,
                                                    IPs.Pressure_Gradient, 
                                                    FIXED_TEMPERATURE_WALL );
                } /* endif */

            } else {
               dX = Grid.xfaceTop(i,j,k)- Grid.Cell[i][j][k].Xc;
               Wl = W[i][j][k] + 
                   ( phi[i][j][k]^ dWdx[i][j][k])*dX.x +
                   ( phi[i][j][k]^ dWdy[i][j][k])*dX.y +
                   ( phi[i][j][k]^ dWdz[i][j][k])*dX.z;
               dX = Grid.xfaceBot(i,j,k+1)- Grid.Cell[i][j][k+1].Xc;
               Wr = W[i][j][k+1] +
                   ( phi[i][j][k+1]^ dWdx[i][j][k+1])*dX.x +
                   ( phi[i][j][k+1]^ dWdy[i][j][k+1])*dX.y +
                   ( phi[i][j][k+1]^ dWdz[i][j][k+1])*dX.z;

            } /* endif */
             
	    // Spacing for Preconditioner 
	    if (IPs.Flow_Type != FLOWTYPE_INVISCID && IPs.Preconditioning) {
  	       delta_n = min(min(TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceE(i,j,k)+Grid.AfaceW(i,j,k))),
	                         TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceN(i,j,k)+Grid.AfaceS(i,j,k)))),
                                 TWO*(Grid.Cell[i][j][k].V/(Grid.AfaceTop(i,j,k)+Grid.AfaceBot(i,j,k))));
	    } /* endif */

            // Evaluate inviscid numerical flux
            switch(IPs.i_Flux_Function) {
              case FLUX_FUNCTION_HLLE :
                Flux =  LES3DTF_pState::FluxHLLE_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                break;
              case FLUX_FUNCTION_ROE :
                Flux =  LES3DTF_pState::FluxRoe_n(Wl, Wr, Grid.nfaceTop(i,j,k));
                break;
              case FLUX_FUNCTION_AUSM_PLUS_UP :
		Flux =  LES3DTF_pState::FluxAUSMplus_up_n(Wl, Wr, Grid.nfaceTop(i,j,k));
	        break;
	    } /* endswitch */

            // Add viscous flux in k direction
            Flux -=  LES3DTF_pState::FluxViscous_n(Wl, Wr,
                                                    W[i][j][k], W[i][j][k+1], 
                                                    dWdx[i][j][k], dWdy[i][j][k], dWdz[i][j][k],
                                                    dWdx[i][j][k+1], dWdy[i][j][k+1], dWdz[i][j][k+1],
                                                    Grid.nfaceTop(i,j,k), Grid.Vot(i,j,k), 
                                                    Grid.delta_ot(i,j,k), Grid.Cell[i][j][k].V, 
                                                    Grid.Cell[i][j][k+1].V, Flow_Type);

            /* Evaluate cell-averaged solution changes. */
             dUdt[i][j][k][k_residual] -= (IPs.CFL_Number*dt[i][j][k])*
                                          Flux*Grid.AfaceTop(i,j,k)/Grid.Cell[i][j][k].V;
             dUdt[i][j][k+1][k_residual] += (IPs.CFL_Number*dt[i][j][k+1])*
                                            Flux*Grid.AfaceBot(i,j,k+1)/Grid.Cell[i][j][k+1].V;

            /* Save top and bottom face boundary flux. */
//           if (j ==  JCl-1) {
//               FluxS[i] = -Flux* Grid.lfaceS(i, j+1);
//           } else if (j ==  JCu) {
//               FluxN[i] = Flux* Grid.lfaceN(i, j);
//           } /* endif */
             
          } /* endfor */

          dUdt[i][j][KCl-1][k_residual] = U_VACUUM;
          dUdt[i][j][KCu+1][k_residual] = U_VACUUM;
      } /* endfor */
   } /* endfor */

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
int Hexa_Block<LES3DTF_pState,LES3DTF_cState>::
Update_Solution_Multistage_Explicit(const int i_stage,
                                    Input_Parameters<LES3DTF_pState,LES3DTF_cState> &IPs){
   
  int i, j, k,  k_residual, index;
  double omega;
  int num_var_update = NumVar()-2;  // Don't update TF and WF using multistage scheme
   

  // Memory for linear system solver.
  //   LES3DTF_cState dU_precon;

  /* Additional variables for dual time stepping. */
  double dTime = ZERO;          // Physical time step
  double residual_denominator;  // Improves convergence for inner iterations
  double length_n = max(max(IPs.Grid_IP.Box_Length, IPs.Grid_IP.Box_Width), IPs.Grid_IP.Box_Height);    

  double local_progress_variable, Yf_u, Yf_b, Yf;
  double lapl_vor, cell_size;

//   if (IPs.Dual_Time_Stepping) {
//     dTime = IPs.dTime;
//   }

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

   for ( k  =  KCl ; k <=  KCu ; ++k ) {
      for ( j  =  JCl ; j <=  JCu ; ++j ) {
         for ( i =  ICl ; i <=  ICu ; ++i ) {
            
  	    // Update conserved solution state
            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING || 
                IPs.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      
	      //U[i][j][k] = Uo[i][j][k] + omega* dUdt[i][j][k][k_residual];

	      for (index=1; index<num_var_update; ++index) {
		U[i][j][k][index] = Uo[i][j][k][index] + omega* dUdt[i][j][k][k_residual][index];
	      }

	      // update last species using: c_n = 1 - sum(1 to N-1) cs
	      U[i][j][k][num_var_update] = U[i][j][k].rho*(ONE - U[i][j][k].sum_species());

	      if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY){
		U[i][j][k].rhok = ZERO;
	      } /* endif */

            } /* endif */

            // Check physical validity of update solution state
            if (IPs.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Mass Fractions, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] << "\n";
		return (1);
              } /* endif */

            } else {
	      if (!U[i][j][k].Realizable_Solution_Check()) {
                cout << "\n " << CFFC_Name() 
                     << " ERROR: Negative Density, Mass Fractions, Kinetic Energy, and/or Sensible Energy: \n"
                     << " cell = (" << i << ", " << j <<", "<< k << ") " 
                     << " X = " <<  Grid.Cell[i][j][k].Xc 
                     << "\n U = " <<  U[i][j][k] 
                     << "\n dUdt = " << dUdt[i][j][k][k_residual] << "\n";
		return (1);
              } /* endif */
	    } /* endif */

            // Update primitive solution state
            W[i][j][k] = U[i][j][k].W();


        
        Yf = W[i][j][k].spec[0].c;
        Yf_u = 0.05518;
        Yf_b = 0.0;
        local_progress_variable = (Yf - Yf_u)/(Yf_b - Yf_u);

	if(local_progress_variable > 0.02  &&  local_progress_variable < 0.98)  {	  

// 	  U[i][j][k].flame.TF_function(U[i][j][k].rhospec[0].c/U[i][j][k].rho,
// 				       U[i][j][k].rhospec[1].c/U[i][j][k].rho,
// 				       U[i][j][k].T(),
// 				       U[i][j][k]._TFactor);

// 	  U[i][j][k].flame.TF_function(U[i][j][k]._TFactor,
// 				       local_progress_variable);

// 	  U[i][j][k].flame.thickening_factor(Grid.Cell[i][j][k].V, U[i][j][k]._laminar_flame_thickness);



	  U[i][j][k].flame.TF = U[i][j][k]._TFactor; // Maximum thickening factor 
  
	  cell_size = pow(Grid.Cell[i][j][k].V, 1.0/3.0);  
 	  lapl_vor = Laplacian_of_Vorticity(*this, i, j, k);
	  U[i][j][k].flame.wrinkling_factor(U[i][j][k]._laminar_flame_speed,
					    U[i][j][k]._laminar_flame_thickness,
					    cell_size,
					    lapl_vor,
					    U[i][j][k].rho,
					    U[i][j][k].mu());

	  if (local_progress_variable > 0.47  &&  local_progress_variable < 0.53)  {
	    U[i][j][k].flame.iso_c_05 = true;
	  }

 	} else {
	  U[i][j][k].flame.WF = ONE;
	  U[i][j][k].flame.TF = U[i][j][k]._TFactor;   
	  U[i][j][k].flame.iso_c_05 = false;

 	} /*endif */

       /* Update power-law variables from the conserved solution */
	W[i][j][k].flame = U[i][j][k].flame;



        /************ FORM LHS FOR DUAL TIME STEPPING SIMI-IMPLICIT WITH PRECONDITIONING ***************/
//         if (Input_Parameters.Local_Time_Stepping == DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
//           double dual_coef_I;
// 	  dual_coef_I = (Input_Parameters.first_step  ?  ONE:THREE/TWO);
          
//           // dSdU
//           dSdU = SemiImplicitBlockJacobi(SolnBlk,i, j,Input_Parameters.Simple_Chemistry);
//           // I
// 	  LinSys.A.identity();
//           // P
//           SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner(Precon,
// 							  SolnBlk.Flow_Type,
// 							  delta_n,
//                                                           length_n,
//                                                           dTime,
//                                                           Input_Parameters.Simple_Chemistry);
   
//           // SolnBlk.Uo[i][j].Low_Mach_Number_Preconditioner_Inverse(Precon_Inv,
// // 	         						 SolnBlk.Flow_Type,
// // 								 delta_n,
// //                                                                  Input_Parameters.Simple_Chemistry);

//           //                             dt
//           //  omega*CFL* dual_coef_I * ------- * I
//           //                            dTime
//           LinSys.A = dual_coef_I*(/*omega*/Input_Parameters.CFL_Number*SolnBlk.dt[i][j]/dTime)*LinSys.A;
//           //                            dt
//           //  omega*CFL*dual_coef_I *  ----- * I - omega*CFL*dt*dSdU 
//           //                           dTime
             
//           LinSys.A -= (/*omega*/Input_Parameters.CFL_Number*SolnBlk.dt[i][j]*dSdU);
    
//           //                              dt
//           //  P + omega*CFL*dual_coef_I  ----- * I - omega*CFL*dt*dSdU 
//           //                             dTime
//           LinSys.A += Precon;   
// 	}
         } /* endfor */    	 
      } /* endfor */    
   } /* endfor */

   /* Solution successfully updated. */
    
   return (0);   

}


/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j,k) of *
 * the computational mesh for the specified             *
 * Hexahedral solution block.  A least squares          *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
template<>
void Hexa_Block<LES3DTF_pState, LES3DTF_cState>::
Linear_Reconstruction_LeastSquares(const int i,
				   const int j,
				   const int k,
				   const int Limiter) {

   int n, n2, n_pts, i_index[26], j_index[26], k_index[26];
   double u0Min, u0Max, uHexa[6], PHI;
   double DxDx_ave, DxDy_ave, DyDy_ave, DxDz_ave, DyDz_ave, DzDz_ave;
   double D;
   Vector3D dX;
   
   /* solnvec in  DU (DUDx_ave, DUDy_ave, DUDz_ave, D1, D2, D3)
      is allocated using new  */
   LES3DTF_pState DU, DUDx_ave, DUDy_ave, DUDz_ave;
   LES3DTF_pState D1, D2, D3;
   LES3DTF_pState Temp1, Temp2, Temp3;
   
   int num_vars = NumVar()-2;  // Don't reconstruct TF and WF
   
   Vector3D dX_neigbor;

   int num=0;
   
   if (i ==  ICl- Nghost || i ==  ICu+ Nghost ||
       j ==  JCl- Nghost || j ==  JCu+ Nghost ||
       k ==  KCl- Nghost || k ==  KCu+ Nghost) {
      n_pts = 0;
   } else {
      n_pts = 26;
      // k plane
      i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
      i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
      i_index[2] = i+1; j_index[2] = j-1; k_index[2] = k;
      i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k;
      i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k;
      i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k;
      i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k;
      i_index[7] = i+1; j_index[7] = j+1; k_index[7] = k;
      //k-1 plane
      i_index[8] = i-1; j_index[8] = j-1; k_index[8] = k-1;
      i_index[9] = i  ; j_index[9] = j-1; k_index[9] = k-1;
      i_index[10] = i+1; j_index[10] = j-1; k_index[10] = k-1;
      i_index[11] = i-1; j_index[11] = j  ; k_index[11] = k-1;
      i_index[12] = i  ; j_index[12] = j  ; k_index[12] = k-1;
      i_index[13] = i+1; j_index[13] = j  ; k_index[13] = k-1;
      i_index[14] = i-1; j_index[14] = j+1; k_index[14] = k-1;
      i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k-1;
      i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k-1;
      //k+1 plane
      i_index[17] = i-1; j_index[17] = j-1; k_index[17] = k+1;
      i_index[18] = i  ; j_index[18] = j-1; k_index[18] = k+1;
      i_index[19] = i+1; j_index[19] = j-1; k_index[19] = k+1;
      i_index[20] = i-1; j_index[20] = j  ; k_index[20] = k+1;
      i_index[21] = i  ; j_index[21] = j  ; k_index[21] = k+1;
      i_index[22] = i+1; j_index[22] = j  ; k_index[22] = k+1;
      i_index[23] = i-1; j_index[23] = j+1; k_index[23] = k+1;
      i_index[24] = i  ; j_index[24] = j+1; k_index[24] = k+1;
      i_index[25] = i+1; j_index[25] = j+1; k_index[25] = k+1;
   } /* endif */
     
   if (n_pts > 0) {
      DUDx_ave.Vacuum();
      DUDy_ave.Vacuum();
      DUDz_ave.Vacuum();
      D1.Vacuum();
      D2.Vacuum();
      D3.Vacuum();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DxDz_ave = ZERO;
      DyDy_ave = ZERO;
      DyDz_ave = ZERO;
      DzDz_ave = ZERO;
      D = ZERO;
      
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
         dX =  Grid.Cell[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ].Xc -
            Grid.Cell[i][j][k].Xc;
         DU =  W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] -  W[i][j][k];
         
         DUDx_ave += DU*dX.x;
         DUDy_ave += DU*dX.y;
         DUDz_ave += DU*dX.z;
         DxDx_ave += dX.x*dX.x;
         DxDy_ave += dX.x*dX.y;
         DxDz_ave += dX.x*dX.z;
         DyDy_ave += dX.y*dX.y;
         DyDz_ave += dX.y*dX.z;
         DzDz_ave += dX.z*dX.z;
         
      } /* endfor */
      
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DUDz_ave = DUDz_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DxDz_ave = DxDz_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      DyDz_ave = DyDz_ave/double(n_pts);
      DzDz_ave = DzDz_ave/double(n_pts);
      
      // (1) Either write a linear solver for 3x3 linear system
      // (2) Or simplely use cramer's rule for this simple system

      D = DxDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
          DxDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D1 = DUDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
           DUDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
           DUDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D2 =DxDx_ave*(DUDy_ave* DzDz_ave - DUDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DUDz_ave - DUDx_ave*DzDz_ave)+
          DxDz_ave*(DUDx_ave*DyDz_ave - DxDz_ave*DUDy_ave);

      D3 =DxDx_ave*(DyDy_ave* DUDz_ave - DyDz_ave*DUDy_ave) +
          DxDy_ave*(DUDx_ave*DyDz_ave - DxDy_ave*DUDz_ave)+
          DxDz_ave*(DxDy_ave*DUDy_ave - DUDx_ave*DyDy_ave);

      dWdx[i][j][k] = D1/D;
      dWdy[i][j][k] = D2/D;
      dWdz[i][j][k] = D3/D;
      
      if (! Freeze_Limiter) {
         for ( n = 1 ; n <= num_vars ; ++n ) {
            
            u0Min =  W[i][j][k][n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               u0Min = min(u0Min,
                           W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ][n]);
               u0Max = max(u0Max,
                           W[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ][n]);
            } /* endfor */
            
            dX =  Grid.xfaceE(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[0] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceW(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[1] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceN(i, j, k)- Grid.Cell[i][j][k].Xc;
            uHexa[2] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
	    dX =  Grid.xfaceS(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[3] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceTop(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[4] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
            dX =  Grid.xfaceBot(i, j, k)- Grid.Cell[i][j][k].Xc;
	    uHexa[5] =  W[i][j][k][n] +
                dWdx[i][j][k][n]*dX.x +
                dWdy[i][j][k][n]*dX.y +
                dWdz[i][j][k][n]*dX.z ;
	    
	    switch(Limiter) {
	    case LIMITER_ONE :
               PHI = ONE;
               break;
	    case LIMITER_ZERO :
               PHI = ZERO;
               break;
	    case LIMITER_BARTH_JESPERSEN :
               PHI = Limiter_BarthJespersen(uHexa,  W[i][j][k][n],
                                            u0Min, u0Max, 6);
               break;
	    case LIMITER_VENKATAKRISHNAN :
               PHI = Limiter_Venkatakrishnan(uHexa,  W[i][j][k][n],
                                             u0Min, u0Max, 6);
               break;
	    case LIMITER_VANLEER :
               PHI = Limiter_VanLeer(uHexa,  W[i][j][k][n],
                                     u0Min, u0Max, 6);
               break;
	    case LIMITER_VANALBADA :
               PHI = Limiter_VanAlbada(uHexa,  W[i][j][k][n],
                                       u0Min, u0Max, 6);
               break;
	    default:
               PHI = Limiter_BarthJespersen(uHexa,  W[i][j][k][n],
                                            u0Min, u0Max, 6);
               break;
 	    }// endswitch
	    
	     phi[i][j][k][n] = PHI;
                        
         } /* endfor */
      } /* endif */
   } else {
       dWdx[i][j][k].Vacuum();
       dWdy[i][j][k].Vacuum();
       dWdz[i][j][k].Vacuum();
       phi[i][j][k].Vacuum();
   } /* endif */
    
}


/********************************************************
 *       Magnitude of the Laplacian of vorticity        *
 ********************************************************/
template<typename HEXA_BLOCK>
double Laplacian_of_Vorticity(HEXA_BLOCK &Solution_Block,
			      const int &i,
			      const int &j,
			      const int &k) {

  Vector3D V, V_neighbour, dVdx_face, dVdy_face, dVdz_face;
  Vector3D dVdx, dVdy, dVdz, dVdx_neighbour, dVdy_neighbour, dVdz_neighbour;

  double Aface, Volume;
  Vector3D Vorticity_Laplacian;
  

  V = Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k],
					  Solution_Block.dWdy[i][j][k],
					  Solution_Block.dWdz[i][j][k]);

  Gradient_of_Vorticity(Solution_Block, i, j, k,
			dVdx, dVdy, dVdz);

  Volume = Solution_Block.Grid.Cell[i][j][k].V;



  //-----------------------------------------------------------//
  //  - Calculate the gradient of vorticity at the cell faces  //
  //                                                           //
  //  - Compute the Laplacian using Green's Theorem            //
  //-----------------------------------------------------------//



  // East face gradient
  V_neighbour = Solution_Block.W[i+1][j][k].vorticity(Solution_Block.dWdx[i+1][j][k],
						      Solution_Block.dWdy[i+1][j][k],
						      Solution_Block.dWdz[i+1][j][k]);

  Gradient_of_Vorticity(Solution_Block, i+1, j, k,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);
  
  Face_Gradient(V, V_neighbour, 
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceE(i,j,k), Solution_Block.Grid.Voe(i,j,k),  
		Solution_Block.Grid.delta_oe(i,j,k), Solution_Block.Grid.Cell[i][j][k].V, 
		Solution_Block.Grid.Cell[i+1][j][k].V,
		dVdx_face, dVdy_face, dVdz_face);


  // East face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceE(i,j,k);

  Vorticity_Laplacian.x = Aface * (dVdx_face.x*Solution_Block.Grid.nfaceE(i,j,k).x +
				   dVdy_face.x*Solution_Block.Grid.nfaceE(i,j,k).y +
				   dVdz_face.x*Solution_Block.Grid.nfaceE(i,j,k).z);

  Vorticity_Laplacian.y = Aface * (dVdx_face.y*Solution_Block.Grid.nfaceE(i,j,k).x +
				   dVdy_face.y*Solution_Block.Grid.nfaceE(i,j,k).y +
				   dVdz_face.y*Solution_Block.Grid.nfaceE(i,j,k).z);

  Vorticity_Laplacian.z = Aface * (dVdx_face.z*Solution_Block.Grid.nfaceE(i,j,k).x +
				   dVdy_face.z*Solution_Block.Grid.nfaceE(i,j,k).y +
				   dVdz_face.z*Solution_Block.Grid.nfaceE(i,j,k).z);



  // West face gradient
  V_neighbour = Solution_Block.W[i-1][j][k].vorticity(Solution_Block.dWdx[i-1][j][k],
						      Solution_Block.dWdy[i-1][j][k],
						      Solution_Block.dWdz[i-1][j][k]);

  Gradient_of_Vorticity(Solution_Block, i-1, j, k,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);

  Face_Gradient(V, V_neighbour,
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceW(i,j,k), Solution_Block.Grid.Vow(i,j,k),
		Solution_Block.Grid.delta_ow(i,j,k), Solution_Block.Grid.Cell[i][j][k].V,
		Solution_Block.Grid.Cell[i-1][j][k].V,
		dVdx_face, dVdy_face, dVdz_face);


  // West face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceW(i,j,k);

  Vorticity_Laplacian.x += Aface * (dVdx_face.x*Solution_Block.Grid.nfaceW(i,j,k).x +
				    dVdy_face.x*Solution_Block.Grid.nfaceW(i,j,k).y +
				    dVdz_face.x*Solution_Block.Grid.nfaceW(i,j,k).z);

  Vorticity_Laplacian.y += Aface * (dVdx_face.y*Solution_Block.Grid.nfaceW(i,j,k).x +
				    dVdy_face.y*Solution_Block.Grid.nfaceW(i,j,k).y +
				    dVdz_face.y*Solution_Block.Grid.nfaceW(i,j,k).z);

  Vorticity_Laplacian.z += Aface * (dVdx_face.z*Solution_Block.Grid.nfaceW(i,j,k).x +
				    dVdy_face.z*Solution_Block.Grid.nfaceW(i,j,k).y +
				    dVdz_face.z*Solution_Block.Grid.nfaceW(i,j,k).z);



  // North face gradient
  V_neighbour = Solution_Block.W[i][j+1][k].vorticity(Solution_Block.dWdx[i][j+1][k],
						      Solution_Block.dWdy[i][j+1][k],
						      Solution_Block.dWdz[i][j+1][k]);

  Gradient_of_Vorticity(Solution_Block, i, j+1, k,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);

  Face_Gradient(V, V_neighbour,
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceN(i,j,k), Solution_Block.Grid.Von(i,j,k),
		Solution_Block.Grid.delta_on(i,j,k), Solution_Block.Grid.Cell[i][j][k].V,
		Solution_Block.Grid.Cell[i][j+1][k].V,
		dVdx_face, dVdy_face, dVdz_face);


  // North face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceN(i,j,k);

  Vorticity_Laplacian.x += Aface * (dVdx_face.x*Solution_Block.Grid.nfaceN(i,j,k).x +
				    dVdy_face.x*Solution_Block.Grid.nfaceN(i,j,k).y +
				    dVdz_face.x*Solution_Block.Grid.nfaceN(i,j,k).z);

  Vorticity_Laplacian.y += Aface * (dVdx_face.y*Solution_Block.Grid.nfaceN(i,j,k).x +
				    dVdy_face.y*Solution_Block.Grid.nfaceN(i,j,k).y +
				    dVdz_face.y*Solution_Block.Grid.nfaceN(i,j,k).z);

  Vorticity_Laplacian.z += Aface * (dVdx_face.z*Solution_Block.Grid.nfaceN(i,j,k).x +
				    dVdy_face.z*Solution_Block.Grid.nfaceN(i,j,k).y +
				    dVdz_face.z*Solution_Block.Grid.nfaceN(i,j,k).z);



  // South face gradient
  V_neighbour = Solution_Block.W[i][j-1][k].vorticity(Solution_Block.dWdx[i][j-1][k],
						      Solution_Block.dWdy[i][j-1][k],
						      Solution_Block.dWdz[i][j-1][k]);

  Gradient_of_Vorticity(Solution_Block, i, j-1, k,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);

  Face_Gradient(V, V_neighbour,
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceS(i,j,k), Solution_Block.Grid.Vos(i,j,k),
		Solution_Block.Grid.delta_os(i,j,k), Solution_Block.Grid.Cell[i][j][k].V,
		Solution_Block.Grid.Cell[i][j-1][k].V,
		dVdx_face, dVdy_face, dVdz_face);


  // South face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceS(i,j,k);

  Vorticity_Laplacian.x += Aface * (dVdx_face.x*Solution_Block.Grid.nfaceS(i,j,k).x +
				    dVdy_face.x*Solution_Block.Grid.nfaceS(i,j,k).y +
				    dVdz_face.x*Solution_Block.Grid.nfaceS(i,j,k).z);

  Vorticity_Laplacian.y += Aface * (dVdx_face.y*Solution_Block.Grid.nfaceS(i,j,k).x +
				    dVdy_face.y*Solution_Block.Grid.nfaceS(i,j,k).y +
				    dVdz_face.y*Solution_Block.Grid.nfaceS(i,j,k).z);

  Vorticity_Laplacian.z += Aface * (dVdx_face.z*Solution_Block.Grid.nfaceS(i,j,k).x +
				    dVdy_face.z*Solution_Block.Grid.nfaceS(i,j,k).y +
				    dVdz_face.z*Solution_Block.Grid.nfaceS(i,j,k).z);



  // Top face gradient
  V_neighbour = Solution_Block.W[i][j][k+1].vorticity(Solution_Block.dWdx[i][j][k+1],
						      Solution_Block.dWdy[i][j][k+1],
						      Solution_Block.dWdz[i][j][k+1]);

  Gradient_of_Vorticity(Solution_Block, i, j, k+1,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);

  Face_Gradient(V, V_neighbour,
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceTop(i,j,k), Solution_Block.Grid.Vot(i,j,k),
		Solution_Block.Grid.delta_ot(i,j,k), Solution_Block.Grid.Cell[i][j][k].V,
		Solution_Block.Grid.Cell[i][j][k+1].V,
		dVdx_face, dVdy_face, dVdz_face);


  // Top face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceTop(i,j,k);

  Vorticity_Laplacian.x += Aface * (dVdx_face.x*Solution_Block.Grid.nfaceTop(i,j,k).x +
				    dVdy_face.x*Solution_Block.Grid.nfaceTop(i,j,k).y +
				    dVdz_face.x*Solution_Block.Grid.nfaceTop(i,j,k).z);

  Vorticity_Laplacian.y += Aface * (dVdx_face.y*Solution_Block.Grid.nfaceTop(i,j,k).x +
				    dVdy_face.y*Solution_Block.Grid.nfaceTop(i,j,k).y +
				    dVdz_face.y*Solution_Block.Grid.nfaceTop(i,j,k).z);

  Vorticity_Laplacian.z += Aface * (dVdx_face.z*Solution_Block.Grid.nfaceTop(i,j,k).x +
				    dVdy_face.z*Solution_Block.Grid.nfaceTop(i,j,k).y +
				    dVdz_face.z*Solution_Block.Grid.nfaceTop(i,j,k).z);



  // Bottom face gradient
  V_neighbour = Solution_Block.W[i][j][k-1].vorticity(Solution_Block.dWdx[i][j][k-1],
						      Solution_Block.dWdy[i][j][k-1],
						      Solution_Block.dWdz[i][j][k-1]);

  Gradient_of_Vorticity(Solution_Block, i, j, k-1,
			dVdx_neighbour, dVdy_neighbour, dVdz_neighbour);

  Face_Gradient(V, V_neighbour,
		dVdx, dVdy, dVdz,
		dVdx_neighbour, dVdy_neighbour, dVdz_neighbour,
		Solution_Block.Grid.nfaceBot(i,j,k), Solution_Block.Grid.Vob(i,j,k),
		Solution_Block.Grid.delta_ob(i,j,k), Solution_Block.Grid.Cell[i][j][k].V,
		Solution_Block.Grid.Cell[i][j][k-1].V,
		dVdx_face, dVdy_face, dVdz_face);


  // Bottom face contribution to the Laplacian
  Aface = Solution_Block.Grid.AfaceBot(i,j,k);

  Vorticity_Laplacian.x += Aface * (dVdx_face.x*Solution_Block.Grid.nfaceBot(i,j,k).x +
				    dVdy_face.x*Solution_Block.Grid.nfaceBot(i,j,k).y +
				    dVdz_face.x*Solution_Block.Grid.nfaceBot(i,j,k).z);

  Vorticity_Laplacian.y += Aface * (dVdx_face.y*Solution_Block.Grid.nfaceBot(i,j,k).x +
				    dVdy_face.y*Solution_Block.Grid.nfaceBot(i,j,k).y +
				    dVdz_face.y*Solution_Block.Grid.nfaceBot(i,j,k).z);

  Vorticity_Laplacian.z += Aface * (dVdx_face.z*Solution_Block.Grid.nfaceBot(i,j,k).x +
				    dVdy_face.z*Solution_Block.Grid.nfaceBot(i,j,k).y +
				    dVdz_face.z*Solution_Block.Grid.nfaceBot(i,j,k).z);


  // Laplacian of vorticity
  Vorticity_Laplacian.x /= Volume;
  Vorticity_Laplacian.y /= Volume; 
  Vorticity_Laplacian.z /= Volume; 


  return  Vorticity_Laplacian.abs();

}


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
			   Vector3D &dVdz) {
  
   int n, n2, n_pts, i_index[26], j_index[26], k_index[26];
   double DxDx_ave, DxDy_ave, DyDy_ave, DxDz_ave, DyDz_ave, DzDz_ave;
   double D;
   
   Vector3D  DU, DUDx_ave, DUDy_ave, DUDz_ave;
   Vector3D  D1, D2, D3;
   Vector3D  dX;


   if (i ==  Solution_Block.ICl- Solution_Block.Nghost || 
       i ==  Solution_Block.ICu + Solution_Block.Nghost ||
       j ==  Solution_Block.JCl- Solution_Block.Nghost || 
       j ==  Solution_Block.JCu + Solution_Block.Nghost ||
       k ==  Solution_Block.KCl - Solution_Block.Nghost || 
       k ==  Solution_Block.KCu + Solution_Block.Nghost) {
     n_pts = 0;
   } else {
     n_pts = 26;
     // k plane
     i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
     i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
     i_index[2] = i+1; j_index[2] = j-1; k_index[2] = k;
     i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k;
     i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k;
     i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k;
     i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k;
     i_index[7] = i+1; j_index[7] = j+1; k_index[7] = k;
     //k-1 plane
     i_index[8] = i-1; j_index[8] = j-1; k_index[8] = k-1;
     i_index[9] = i  ; j_index[9] = j-1; k_index[9] = k-1;
     i_index[10] = i+1; j_index[10] = j-1; k_index[10] = k-1;
     i_index[11] = i-1; j_index[11] = j  ; k_index[11] = k-1;
     i_index[12] = i  ; j_index[12] = j  ; k_index[12] = k-1;
     i_index[13] = i+1; j_index[13] = j  ; k_index[13] = k-1;
     i_index[14] = i-1; j_index[14] = j+1; k_index[14] = k-1;
     i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k-1;
     i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k-1;
     //k+1 plane
     i_index[17] = i-1; j_index[17] = j-1; k_index[17] = k+1;
     i_index[18] = i  ; j_index[18] = j-1; k_index[18] = k+1;
     i_index[19] = i+1; j_index[19] = j-1; k_index[19] = k+1;
     i_index[20] = i-1; j_index[20] = j  ; k_index[20] = k+1;
     i_index[21] = i  ; j_index[21] = j  ; k_index[21] = k+1;
     i_index[22] = i+1; j_index[22] = j  ; k_index[22] = k+1;
     i_index[23] = i-1; j_index[23] = j+1; k_index[23] = k+1;
     i_index[24] = i  ; j_index[24] = j+1; k_index[24] = k+1;
     i_index[25] = i+1; j_index[25] = j+1; k_index[25] = k+1;
   }

   if (n_pts > 0) {
      DUDx_ave.zero();
      DUDy_ave.zero();
      DUDz_ave.zero();
      D1.zero();
      D2.zero();
      D3.zero();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DxDz_ave = ZERO;
      DyDy_ave = ZERO;
      DyDz_ave = ZERO;
      DzDz_ave = ZERO;
      D = ZERO;
      
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
         dX = Solution_Block.Grid.Cell[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ].Xc 
	      - Solution_Block.Grid.Cell[i][j][k].Xc;

	 DU = Solution_Block.W[ i_index[n2] ]
	                      [ j_index[n2] ]
                              [ k_index[n2] ].vorticity(Solution_Block.dWdx[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ],
							Solution_Block.dWdy[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ],
							Solution_Block.dWdz[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ]) 
	      - Solution_Block.W[i][j][k].vorticity(Solution_Block.dWdx[i][j][k],
						    Solution_Block.dWdy[i][j][k],
						    Solution_Block.dWdz[i][j][k]);

         DUDx_ave += DU*dX.x;
         DUDy_ave += DU*dX.y;
         DUDz_ave += DU*dX.z;
         DxDx_ave += dX.x*dX.x;
         DxDy_ave += dX.x*dX.y;
         DxDz_ave += dX.x*dX.z;
         DyDy_ave += dX.y*dX.y;
         DyDz_ave += dX.y*dX.z;
         DzDz_ave += dX.z*dX.z;
         
      } /* endfor */
      
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DUDz_ave = DUDz_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DxDz_ave = DxDz_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      DyDz_ave = DyDz_ave/double(n_pts);
      DzDz_ave = DzDz_ave/double(n_pts);
     
   
      // use cramer's rule for this simple system

      D = DxDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
          DxDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D1 = DUDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
           DUDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
           DUDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D2 =DxDx_ave*(DUDy_ave* DzDz_ave - DUDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DUDz_ave - DUDx_ave*DzDz_ave)+
          DxDz_ave*(DUDx_ave*DyDz_ave - DxDz_ave*DUDy_ave);

      D3 =DxDx_ave*(DyDy_ave* DUDz_ave - DyDz_ave*DUDy_ave) +
          DxDy_ave*(DUDx_ave*DyDz_ave - DxDy_ave*DUDz_ave)+
          DxDz_ave*(DxDy_ave*DUDy_ave - DUDx_ave*DyDy_ave);

      dVdx = D1/D;
      dVdy = D2/D;
      dVdz = D3/D;  
   } else {
      dVdx.zero();
      dVdy.zero();
      dVdz.zero();
   } /* endif */
   
}


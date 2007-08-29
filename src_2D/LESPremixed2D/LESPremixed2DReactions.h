/****************** LESPremixedReactions.h **************************************
  This class defines the specializations of the Reaction class.

  
          - based on Chem2D -> Reactions.h

***********************************************************************/
#ifndef _LESPREMIXED2D_REACTIONS_INCLUDED
#define _LESPREMIXED2D_REACTIONS_INCLUDED



/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using the general law of mass action.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations
              of LESPremixed2D_pState. (c_i*rho/M_i)   mol/m^3

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
template<>
inline void Reaction_set::omega<LESPremixed2D_pState,LESPremixed2D_cState>
                               (LESPremixed2D_cState &U,
			        const LESPremixed2D_pState &W,  
			        const int Flow_Type ) const {
 
  double Temp = W.T();  // K
  double rho= W.rho/THOUSAND; // kg/m^3 -> g/cm^3
  double pressure = W.p*TEN;  // N/m^2 -> dyne/cm^2
  double a,b, ans(ZERO);

  double *Wdot = NULL;
  if (reactset_flag == CH4_15STEP_ARM2  ||  reactset_flag == CH4_15STEP_ARM3) {
    Wdot = new double[num_species];
  }

  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND;  //kg/mol -> g/mol
    c[i] = W.spec[i].c;                        //unitless
  }

  switch(reactset_flag){
  case NO_REACTIONS:
    cerr<<"\n You shouldn't get here, NO_REACTIONS in Reactionset::omeag(..)";
    exit(1);
    break;
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
  case CH4_1STEP: 
    //laminar case     
    a=1.0;//   a = 0.2; 
    b=1.0;//     b = 1.3;
    kf[0] = reactions[0].kf(Temp)*pow((W.SpecCon(0))/MILLION,a)*pow((W.SpecCon(1)/MILLION),b);        
    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //CH4
	ans = - kf[0];
	break;
      case 1 : //O2
	ans = - TWO*kf[0];
	break;
      case 2 : //CO2
	ans = kf[0];
	break;
      case 3 : //H2O
	ans = TWO*kf[0];
	break;
      };
      //ans in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3      
      U.rhospec[index].c = M[index]*ans*THOUSAND;  
    }
  
    break;
    
    //TWO STEP CH4
  case CH4_2STEP:
    kf[0] = reactions[0].kf(Temp)*pow(W.SpecCon(0)/MILLION,0.2)*pow(W.SpecCon(1)/MILLION,1.3);
    kb[0] = ZERO;
    kf[1] = reactions[1].kf(Temp)*(W.SpecCon(4)/MILLION)*pow(W.SpecCon(3)/MILLION,0.5)*pow(W.SpecCon(1)/MILLION,0.25); 
    kb[1] = reactions[1].kb(Temp)*(W.SpecCon(2)/MILLION);
   
    //cout<<"\n Sw "<<Temp<<" "<<reactions[0].kf(Temp)<<" "<<reactions[1].kf(Temp)<<" "<<reactions[1].kb(Temp);

    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //CH4
	ans = - kf[0];
	break;
      case 1:  //O2
	ans = -1.5*kf[0] - HALF*kf[1] + HALF*kb[1];
	break;
      case 2:  //CO2
	ans = kf[1]  - kb[1];
	break;
      case 3:  //H2O
	ans = TWO*kf[0];
	break;
      case 4:  //CO
	ans = kf[0] - kf[1] + kb[1];
 	break;
      };     
      //ans in kg/m^3*s     g/mol *(mol/cm^3*s)*1e3       
      U.rhospec[index].c = M[index]*ans*THOUSAND;
    }
    break;
  
  case H2O2_1STEP:
   
      kf[0] = reactions[0].kf(Temp);
      
      for(int index =0; index<num_react_species; index++){
	switch(index) {
	case 0 : //H2
	  ans = - TWO*kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	case 1 : //O2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	case 2 : //H2O
	  ans = TWO*kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	};
	//ans in kg/m^3*s   
	U.rhospec[index].c = M[index]*ans*THOUSAND;
      }
      break;
  
      //TWO STEP H2&O2
  case H2O2_2STEP:
   
      kf[0] = reactions[0].kf(Temp,W.spec[0].c,W.spec[1].c,W.spec[4].c);
      kb[0] = kf[0]/(reactions[0].keq(W,Temp));
      kf[1] = reactions[1].kf(Temp,W.spec[0].c,W.spec[1].c,W.spec[4].c);
      kb[1] = kf[1]/(reactions[1].keq(W,Temp));
      //cout<<"\n kf1 "<< kf[0]<<" kb1 "<< kb[0]<<" kf2 "<<kf[1]<<" kb2 "<< kb[1];
      
      for(int index =0; index<num_react_species; index++){
	switch(index) {
	case 0 : //H2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 + kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12
	    - kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 + kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12;
	  break;
	case 1 : //O2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 + kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12;
	  break;
	case 2 : //OH
	  ans = TWO*(kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 - kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12)
	    + TWO*( -kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 + kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12);
	  break;
	case 3 : //H2O
	  ans = TWO*(kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 - kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12);      
	  break;
	};
	//ans in kg/m^3*s   
	U.rhospec[index].c = M[index]*ans*THOUSAND;
      }
      break;

    // 8 STEP H2 & O2
  case H2O2_8STEP: 
    
    for(int spec =0; spec<num_reactions; spec++){
      kf[spec] = reactions[spec].kf(Temp);
      kb[spec] = kf[spec]/reactions[spec].keq(W,Temp);
    }
 
    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //O
	ans = (kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]-kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]-kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]+kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]+kf[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])-kb[3]*rho*rho*c[4]/M[4]*c[0]/M[0]-kf[6]
			*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]+kb[6]*rho*rho*c[5]/M[5]*c[4]/M[4]
			-2.0*kf[7]*rho*rho*rho*c[0]*c[0]/(M[0]*M[0])*c[4]/M[4]+2.0*kb[7]*rho*rho*c[1]/M
			[1]*c[4]/M[4]);
	
	break;
      case 1 : //O2
	ans = (-kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]+kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[7]*rho*rho*rho*c[0]*c[0]/(M[0]*M[0])*c[4]/M[4]-kb[7]*rho*rho*c
			[1]/M[1]*c[4]/M[4]);
  
	break;
      case 2 : //H
	ans = (-kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]+kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]-kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]+kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]-kb[2]*rho*rho*c[2]/M[2]*c[4]/M[4]-2.0*kf
			[4]*rho*rho*rho*c[2]*c[2]/(M[2]*M[2])*c[4]/M[4]+2.0*kb[4]*rho*rho*c[3]/M[3]*c
			[4]/M[4]-kf[5]*rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]+kb[5]*rho*rho*c[4]*c
			 [4]/(M[4]*M[4])-kf[6]*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]+kb[6]*rho*rho*c
			[5]/M[5]*c[4]/M[4]);	
	break;
      case 3 : //H2
	ans = (-kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]+kb[1]*rho*rho*c[5]/M
			 [5]*c[2]/M[2]-kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]+kb[2]*rho*rho*c[2]/M[2]*c[4]/M
			 [4]+kf[4]*rho*rho*rho*c[2]*c[2]/(M[2]*M[2])*c[4]/M[4]-kb[4]*rho*rho*c[3]/M[3]*c
			 [4]/M[4]);
	break;
	
      case 4 : //H2O
	ans = (kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]-kb[2]*rho*rho*c[2]/M
			[2]*c[4]/M[4]+kf[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])-kb[3]*rho*rho*c[4]/M[4]*c[0]/
			M[0]+kf[5]*rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]-kb[5]*rho*rho*c[4]*c[4]/(M[4]*M[4]));
	break;

      case 5 : // OH
	ans = (kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]-kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]-kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]-kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]+kb[2]*rho*rho*c[2]/M[2]*c[4]/M[4]-2.0*kf
			[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])+2.0*kb[3]*rho*rho*c[4]/M[4]*c[0]/M[0]-kf[5]*
			rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]+kb[5]*rho*rho*c[4]*c[4]/(M[4]*M[4])+
			kf[6]*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]-kb[6]*rho*rho*c[5]/M[5]*c[4]/M
			[4]);
	break;
      };
      U.rhospec[index].c = M[index]*ans*THOUSAND;
    }
    break;


    // 15 step CH4 based on GRI 2.11
  case CH4_15STEP_ARM2:
    // compute reaction rates calling the subroutine CKWYP15STEP211
    // Wdot - mol/(cm^3*s)
    ckwyp15step211_(pressure, Temp, c, Wdot);

    for (int i=0; i<num_react_species; ++i) {
      // in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3             
      U.rhospec[i].c = M[i]*Wdot[i]*THOUSAND;
    }
    break;


    // 15 step CH4 based on GRI 3
  case CH4_15STEP_ARM3:
    // compute reaction rates calling the subroutine CKWYP15STEP3
    // Wdot - mol/(cm^3*s)
    ckwyp15step30_(pressure, Temp, c, Wdot);

    for (int i=0; i<num_react_species; ++i) {
      // in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3             
      U.rhospec[i].c = M[i]*Wdot[i]*THOUSAND;
    }
    break;

    //---------------------------------//
    //----- User Specified ------------//
    //---------------------------------//
  case USER:
    cerr<<"\nUser specified not currently available in shareware version";
    exit(1);
    break;
  default:
    break;
  };

  //clean up memory
  if (reactset_flag == CH4_15STEP_ARM2  ||  reactset_flag == CH4_15STEP_ARM3) { 
    delete[] Wdot; 
  }

} //end omega()


/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables by using finite difference
  
   dSwdU:  Matrix of source terms 
************************************************************************/
template<>
inline void Reaction_set::Finite_Difference_dSwdU<LESPremixed2D_pState>
                          (DenseMatrix &dSwdU, 
			   const LESPremixed2D_pState &Wlocal,
			   const int Flow_Type) const {
  
  const double b = numeric_limits<double>::epsilon();  // 1.0E-12;  
  LESPremixed2D_pState W1;
  LESPremixed2D_cState S1, S, Ulocal, U1, EPS;
  
  int lower_index = Wlocal.NUM_VAR_LESPREMIXED2D + 1 - (Wlocal.ns + Wlocal.nscal);
  
  Ulocal = U(Wlocal);
  S = Wlocal.Sw(Wlocal.React.reactset_flag, Flow_Type);

  // Perturbation parameter (EPSILON)
  for (int j=lower_index; j < Wlocal.NUM_VAR_LESPREMIXED2D; ++j) {
    EPS[j] = sqrt( b*fabs(Ulocal[j]) + b )/1000.0;
  }
  
  U1 = Ulocal;
  S1 = S;
  for (int j=lower_index; j < Wlocal.NUM_VAR_LESPREMIXED2D; ++j) {
    if ( EPS[j] != 0.0 ) {
      U1[j] += EPS[j];
      W1 = W(U1);
      S1 = W1.Sw(W1.React.reactset_flag, Flow_Type);    
      for (int i=lower_index; i < Wlocal.NUM_VAR_LESPREMIXED2D; ++i) {
	dSwdU(i-1, j-1) = TwoPointFiniteDifference(S[i], S1[i], EPS[j]);
      }
      // reset U1 to U
      U1 = Ulocal;
    }
  } // end for

} // end  Finite_Difference_dSwdU



/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables by using the complex step method
  
   dSwdU:  Matrix of source terms 
************************************************************************/
template<>
inline void Reaction_set::Complex_Step_dSwdU<LESPremixed2D_pState>
                          (DenseMatrix &dSwdU, 
			   const LESPremixed2D_pState &Wlocal) const {

  const double TOL = 1.0E-15, eps = 1.0E-100;
  cplx *Y, *M, *Wdot_perturbed;

  //assert(num_species == Wlocal.ns);

  // allocate 
  Y = new cplx[num_species];  
  M = new cplx[num_species];
  Wdot_perturbed = new cplx[num_species]; 


  // pressure in dyne/cm^2
  cplx pressure(TEN*Wlocal.p, ZERO);
  // temperature in K 
  cplx temperature(Wlocal.T(), ZERO);

  for (int i=0; i<Wlocal.ns; ++i) {
    Y[i] = Wlocal.spec[i].c;    //( Wlocal.spec[i].c < TOL )  ?  TOL : Wlocal.spec[i].c;
    M[i] = Wlocal.specdata[i].Mol_mass();  // kg/kg-mol
    //Wdot_perturbed[i] = 0.0;
  }

  
  int lower_index = Wlocal.NUM_VAR_LESPREMIXED2D + 1 - Wlocal.ns;
  cplx c_eps(0.0, eps);

  for (int j=lower_index; j < Wlocal.NUM_VAR_LESPREMIXED2D; ++j) {
    if ( Wlocal.spec[j].c > TOL ) {

      Y[j-lower_index] += c_eps;    

      // compute the perturbed reaction rates moles/(cm**3*sec)
      switch(reactset_flag) {
      case CH4_15STEP_ARM2 :
      	cplx15step211_(pressure, temperature, Y, Wdot_perturbed);
	break;
      case CH4_15STEP_ARM3 :
	cplx15step30_(pressure, temperature, Y, Wdot_perturbed);
	break;
      default :
        cplx15step30_(pressure, temperature, Y, Wdot_perturbed);
	break;
      }
  
      for (int i=lower_index; i < Wlocal.NUM_VAR_LESPREMIXED2D; ++i) {
	// in kg/m^3*s  ->  (kg/mol)*(mol/cm^3*s)*1e6 
	Wdot_perturbed[i-lower_index] *= M[i-lower_index]*MILLION;

	dSwdU(i-1, j-1) = imag( Wdot_perturbed[i-lower_index]/eps )/Wlocal.rho;
	//Wdot_perturbed[i-lower_index] = (0.0, 0.0);
      }
      
      // reset Y to its unperturbed value
      Y[j-lower_index] -= c_eps;

    } // end if
  } // end for

  // deallocate
  delete[]  Y;
  delete[]  M;
  delete[]  Wdot_perturbed;

}  // end Complex_Step_dSwdU


/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables.  Ie it returns conserved data.
  
   dSwdU:  Matrix of source terms 
   W:      Primitive State data.

  These are added to the RHS so += is used.

  NOTE: Be carful of mass fractions in the denominator as 
        if they go to ZERO, it will cause floating point exceptions
       (ie. Divison by zero);  
************************************************************************/
template<>
inline void Reaction_set::dSwdU<LESPremixed2D_pState,LESPremixed2D_cState>
                               (DenseMatrix &dSwdU,
				const LESPremixed2D_pState &W, 
				const bool &CFL_flag, 
				const int Flow_Type, 
				const int Solver_type) const {

  /***************** Local Variables *****************/
  double Temp = W.T();
  double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3


//   double Con0, Con1, Rt;
//   double dkf0_dp, dkf0_drho;
//   double dkf0_dc1, dkf0_dc2, dkf0_dc3, dkf0_dc4; 
//   //BEING LAZY INSTEAD OF ANALYTICALY DETERMINING dSdU
//   DenseMatrix dSwdW(W.NUM_VAR_LESPREMIXED2D-1,W.NUM_VAR_LESPREMIXED2D-1,ZERO);  //SHOULD BE MOVED TO TEMP AS WELL!!!
//   DenseMatrix dWdQ(W.NUM_VAR_LESPREMIXED2D-1,W.NUM_VAR_LESPREMIXED2D-1,ZERO);
  
  double VALUE = TOLER; //sqrt(TOLER);
  //////////////////////////////////////////////////////////////////////////////
  // THIS HACK DOESN'T REALLY WORK, ESPECIALLY FOR 2STEP EQUATIONS !!!!!!!!!!!!!
  //////////////////////////////////////////////////////////////////////////////
  if(Solver_type == IMPLICIT){ VALUE=TOLER*TOLER; }

  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND; //kg/mol -> g/mol 
    c[i] = W.spec[i].c;
    
    //For handling ~= ZERO mass fractions that appear in the denominator of dSwdU
    //by setting a lower tolerance allowed, and if below that set to tolerance 
    if( c[i] < VALUE){
      c_denom[i] = VALUE;
    } else {
      c_denom[i] = c[i];
    }
  }

  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W.nscal;
  /*******************************************
   *  Reaction Mechanism Jacobians           *
   *                                         * 
   *******************************************/
  switch(reactset_flag){
  case NO_REACTIONS:
    //this case shouldn't be called
    cerr<<" Ummm, this dSwdU shouldn't be called for with NO_REACTIONS, check quad_singleblock.";
    break;
    
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
    // This is far from an elegant solution, but its easily created using maple.
    // It can probably be simplified for faster computation.

  case CH4_1STEP:
    kf[0] = reactions[0].kf(Temp); 

    // One forward reaction ... so only depending on the fuel and oxidize ..  
    // for coef a=0.2, b=1.3
 //    dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1);
//     dSwdU(NUM_VAR+1,NUM_VAR) += -0.4*M[1]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
//     dSwdU(NUM_VAR+2,NUM_VAR) += 0.2*M[2]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
//     dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];    
//     dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
//     dSwdU(NUM_VAR+1,NUM_VAR+1) += -0.26E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3);
//     dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.13E1*M[2]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
//     dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];  
 
    //for coef a=b=1.0
    dSwdU(NUM_VAR,NUM_VAR) += -1.0*kf[0]*rho*c[1]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+1) += -1.0*kf[0]*rho*c[0]/M[1];
    dSwdU(NUM_VAR+1,NUM_VAR) += -2.0*kf[0]*rho*c[1]/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -2.0*kf[0]*rho*c[0]/M[0];    
    dSwdU(NUM_VAR+2,NUM_VAR) +=  M[2]*kf[0]*rho*c[1]/M[1]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += M[2]*kf[0]*rho*c[0]/M[0]/M[1];
    dSwdU(NUM_VAR+3,NUM_VAR) += 2.0*M[3]*kf[0]*rho*c[1]/M[1]/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+1) += 2.0*M[3]*kf[0]*rho*c[0]/M[0]/M[1];


//     /********** dSwdW ***********************/
//     Con0 = rho*c[0]/M[0];
//     Con1 = rho*c[1]/M[0];
    
//     Rt = W.Rtot(); 
//     dkf0_dp = reactions[0].dkf_dT(Temp)/(rho*Rt); ///THOUSAND;
//     dkf0_drho = reactions[0].dkf_dT(Temp)*(-Temp/rho);
    
//     //take into account kf(rho,p,and ci)
//     dkf0_dc1 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[0].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc2 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[1].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc3 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[2].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc4 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[3].Rs()-W.specdata[4].Rs()))*Con0*Con1;

//     dSwdW(NUM_VAR,NUM_VAR)  += -M[0]*dkf0_dc1;
//     dSwdW(NUM_VAR+1,NUM_VAR) += -TWO*M[1]*dkf0_dc1;
//     dSwdW(NUM_VAR+2,NUM_VAR) += M[2]*dkf0_dc1;
//     dSwdW(NUM_VAR+3,NUM_VAR) += TWO*M[3]*dkf0_dc1;

//     dSwdW(NUM_VAR,NUM_VAR+1)  += -M[0]*dkf0_dc2;
//     dSwdW(NUM_VAR+1,NUM_VAR+1) += -TWO*M[1]*dkf0_dc2;
//     dSwdW(NUM_VAR+2,NUM_VAR+1) += M[2]*dkf0_dc2;
//     dSwdW(NUM_VAR+3,NUM_VAR+1) += TWO*M[3]*dkf0_dc2;
    
//     dSwdW(NUM_VAR,NUM_VAR+2)  += -M[0]*dkf0_dc3;
//     dSwdW(NUM_VAR+1,NUM_VAR+2) += -TWO*M[1]*dkf0_dc3;
//     dSwdW(NUM_VAR+2,NUM_VAR+2) += M[2]*dkf0_dc3;
//     dSwdW(NUM_VAR+3,NUM_VAR+2) += TWO*M[3]*dkf0_dc3;

//     dSwdW(NUM_VAR,NUM_VAR+3)  += -M[0]*dkf0_dc4;
//     dSwdW(NUM_VAR+1,NUM_VAR+3) += -TWO*M[1]*dkf0_dc4;
//     dSwdW(NUM_VAR+2,NUM_VAR+3) += M[2]*dkf0_dc4;
//     dSwdW(NUM_VAR+3,NUM_VAR+3) += TWO*M[3]*dkf0_dc4;

//     dSwdW(NUM_VAR,0) += -M[0]*dkf0_drho*Con0*Con1 - kf[0]*Con1*c[0] - M[0]*kf[0]*Con0*c[1]/M[1];    
//     dSwdW(NUM_VAR+1,0) += -TWO*M[1]*dkf0_drho*Con0*Con1 - TWO*kf[0]*M[1]*Con1*c[0]/M[0] - TWO*kf[0]*Con0*c[1];    
//     dSwdW(NUM_VAR+2,0) += M[2]*dkf0_drho*Con0*Con1 + kf[0]*M[2]*Con1*c[0]/M[0] + M[2]*kf[0]*Con0*c[1]/M[1];    
//     dSwdW(NUM_VAR+3,0) += TWO*M[3]*dkf0_drho*Con0*Con1 + TWO*kf[0]*M[3]*Con1*c[0]/M[0] + TWO*M[3]*kf[0]*Con0*c[1]/M[1];   

//     dSwdW(NUM_VAR,3) += -M[0]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+1,3) += -TWO*M[1]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+2,3) += M[2]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+3,3) += TWO*M[3]*dkf0_dp*Con0*Con1;
   
//     dSwdW(NUM_VAR,NUM_VAR)   += -kf[0]*Con1*rho*THOUSAND;
//     dSwdW(NUM_VAR+1,NUM_VAR) += -TWO*M[1]*kf[0]*Con1*rho/M[0]*THOUSAND;
//     dSwdW(NUM_VAR+2,NUM_VAR) += M[2]*kf[0]*Con1*rho/M[0]*THOUSAND;
//     dSwdW(NUM_VAR+3,NUM_VAR) += TWO*M[3]*kf[0]*Con1*rho/M[0]*THOUSAND;

//     dSwdW(NUM_VAR,NUM_VAR+1)   += -M[0]*kf[0]*Con0*rho/M[1]*THOUSAND;
//     dSwdW(NUM_VAR+1,NUM_VAR+1) += -TWO*kf[0]*Con0*rho*THOUSAND;
//     dSwdW(NUM_VAR+2,NUM_VAR+1) += M[2]*kf[0]*Con0*rho/M[1]*THOUSAND;
//     dSwdW(NUM_VAR+3,NUM_VAR+1) += TWO*M[3]*kf[0]*Con0*rho/M[1]*THOUSAND;

//     // LAZY dSwdU = dSwdW * dWdU
//     W.dWdU(dWdQ,Flow_Type); 
//     dSwdU += dSwdW*dWdQ;
//     /***************************************/

    //this is a work around for the delta t calculation using an unesseccarily small value
    //when c[0] -> ZERO
    if(c_denom[0] != c[0] && CFL_flag){ dSwdU(NUM_VAR,NUM_VAR)=ZERO; }
    
    break;


    //TWO STEP CH4   
  case CH4_2STEP:  
    /******************** ORIGINAL ****************************************/
    //still some issues with units ??? ie.  dSwdU(6,6) ??
    //which is also on the diagonal so messes with CFL???
    kf[0] = reactions[0].kf(Temp);      
    kb[0] = ZERO;
    kf[1] = reactions[1].kf(Temp);  
    kb[1] = reactions[1].kb(Temp);

    //cout<<"\n dSwdU "<<Temp<<" "<<reactions[0].kf(Temp)<<" "<<reactions[1].kf(Temp)<<" "<<reactions[1].kb(Temp);

    dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3);
    dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];      
    
    dSwdU(NUM_VAR+1,NUM_VAR) += -0.3*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -rho*(0.195E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
			*M[4]+0.125*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/M[1]/M[4]/pow(rho*c_denom[1]/M[1],0.75);
    dSwdU(NUM_VAR+1,NUM_VAR+2) += 0.5*M[1]*kb[1]/M[2];
    dSwdU(NUM_VAR+1,NUM_VAR+3) += -0.25*M[1]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+1,NUM_VAR+4) += -0.5*M[1]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
   
    dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.25*M[2]*kf[1]*rho*c[4]/M[4]*sqrt(rho*c[3]/M[3])/pow(rho*c_denom[1]/M[1],0.75)/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -kb[1];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += 0.5*M[2]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+2,NUM_VAR+4) += M[2]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
    
    dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];    
   
    dSwdU(NUM_VAR+4,NUM_VAR) += 0.2*M[4]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+4,NUM_VAR+1) += rho*(0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
		       *M[4]-0.25*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/(M[1]*M[1])/pow(rho*c_denom[1]/M[1],0.75);
    dSwdU(NUM_VAR+4,NUM_VAR+2) += M[4]*kb[1]/M[2];
    dSwdU(NUM_VAR+4,NUM_VAR+3) += -0.5*kf[1]*rho*c[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+4,NUM_VAR+4) += -1.0*kf[1]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);

//     if(!CFL_flag) cout<<"\n ORIGINAL\n"<<dSwdU;
//     dSwdU.zero();

//     /******************** MODIFIED ****************************************/   
//     //with kf(T), kb(T)

//     kf[0] = reactions[0].kf(Temp);      //UNITS ISSUES due to rho and M and kf's
//     kb[0] = ZERO;
//     kf[1] = reactions[1].kf(Temp);  
//     kb[1] = reactions[1].kb(Temp);
    
//     double dkf_dT[2], dkb_dT[2];
//     dkf_dT[0] = reactions[0].dkf_dT(Temp);
//     dkb_dT[0] = ZERO;
//     dkf_dT[1] = reactions[1].dkf_dT(Temp);
//     dkb_dT[1] = reactions[1].dkb_dT(Temp);

//     double T = Temp;
//     double Cp = W.Cp();      
//     double u = W.v.x*HUNDRED; 
//     double v = W.v.y*HUNDRED; 
//     //double p = W.p;      
//     double Rtot = W.Rtot();   
//     double htot = W.h();      
//     double h[6],Rs[6];

//     for(int i=0; i<num_species; i++){    
//       h[i] = (W.specdata[i].Enthalpy(Temp)+W.specdata[i].Heatofform());  
//       Rs[i] = (W.specdata[i].Rs()); 
//     }    

//       if(!CFL_flag) cout<<"\n NEW\n"<<dSwdU;
   
    //this is a work around for the delta t calculation using an unesseccarily small value
    //when CH4 & O2 -> ZERO
    //if( (c[0] < TOLER && CFL_flag) || (c[0] < TOLER && CFL_flag)  ){  
    
    if( c_denom[0] != c[0] && CFL_flag ){
      for(int i=0; i<num_react_species; i++){
	dSwdU(NUM_VAR+i,NUM_VAR+i)=ZERO; 
      }
    }
    break;
  

  case H2O2_1STEP:  
    kf[0] = reactions[0].kf(Temp);
    rho = rho;
    M[0] = M[0];
    M[1] = M[1];

    dSwdU(NUM_VAR,NUM_VAR) += -4.0*rho*rho/M[0]*kf[0]*c[0]*c[1]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+1) = -2.0*rho*rho/M[0]*kf[0]*c[0]*c[0]/M[1];
    
    dSwdU(NUM_VAR+1,NUM_VAR) = -2.0*rho*rho*kf[0]*c[0]/(M[0]*M[0])*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+1) = -rho*rho*kf[0]*c[0]*c[0]/(M[0]*M[0]);
    
    dSwdU(NUM_VAR+2,NUM_VAR) = 4.0*rho*rho*M[2]*kf[0]*c[0]/(M[0]*M[0])*c[1]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+1)= 2.0*rho*rho*M[2]*kf[0]*c[0]*c[0]/(M[0]*M[0])/M[1];
  
    break;
   //TWO STEP H2 and O2
  case H2O2_2STEP:
    kf[0] = reactions[0].kf(Temp,c[0],c[1],W.spec[4].c);
    kb[0] = kf[0]/reactions[0].keq(W,Temp);
    kf[1] = reactions[1].kf(Temp,c[0],c[1],W.spec[4].c);
    kb[1] = kf[1]/reactions[1].keq(W,Temp);  
  
    dSwdU(NUM_VAR,NUM_VAR) += -rho*(kf[0]*c[1]*M[2]*M[2]+kf[1]*rho*c[2]*c[2]*M[1])/M[1]/(M[2]*M[2]);
    dSwdU(NUM_VAR,NUM_VAR+1) += -kf[0]*rho*c[0]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+2) += -2.0*rho*c[2]*(-kb[0]*M[0]+kf[1]*rho*c[0])/(M[2]*M[2]);
    dSwdU(NUM_VAR,NUM_VAR+3) += 2.0*M[0]*kb[1]*rho*c[3]/(M[3]*M[3]);    
  
    dSwdU(NUM_VAR+1,NUM_VAR) += -kf[0]*rho/M[0]*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -kf[0]*rho*c[0]/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+2) += 2.0*M[1]*kb[0]*rho*c[2]/(M[2]*M[2]);    
   
    dSwdU(NUM_VAR+2,NUM_VAR) += rho*(0.2E1*kf[0]*c[1]*M[2]*M[2]-0.2E1*kf[1]*rho*c[2]*c[2]*M[1])/M[2]/M[0]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += 2.0*M[2]*kf[0]*rho*c[0]/M[0]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -rho*c[2]*(0.4E1*kb[0]*M[0]+0.4E1*kf[1]*rho*c[0])/M[2]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += 4.0*M[2]*kb[1]*rho*c[3]/(M[3]*M[3]);    
   
    dSwdU(NUM_VAR+3,NUM_VAR) += 2.0*M[3]*kf[1]*rho*rho*c[2]*c[2]/(M[2]*M[2])/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+2) += 4.0*M[3]*kf[1]*rho*rho*c[2]/(M[2]*M[2])*c[0]/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+3) += -4.0/M[3]*kb[1]*rho*c[3];  

    break;
    // 8 STEP H2 & O2
  case H2O2_8STEP:  

    for(int spec =0; spec<num_reactions; spec++){
      kf[spec] = reactions[spec].kf(Temp);
      kb[spec] = kf[spec]/reactions[spec].keq(W,Temp); 
    }

    dSwdU(NUM_VAR,NUM_VAR) += -1/M[0]*rho*(kb[0]*c[5]*M[0]*M[3]*M[4]*M[2]+kf[1]*c[3]*M[5]
					   *M[0]*M[4]*M[2]+kb[3]*c[4]*M[5]*M[0]*M[3]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[0]*M
			       [3]+4.0*kf[7]*rho*c[0]*c[4]*M[5]*M[3]*M[2])/M[5]/M[3]/M[4]/M[2];
    dSwdU(NUM_VAR,NUM_VAR+1) += M[0]*rho*(kf[0]*c[2]*M[4]+2.0*kb[7]*c[4]*M[2])/M[2]/M[1]/M[4];
    dSwdU(NUM_VAR,NUM_VAR+2) += -rho*(-kf[0]*c[1]*M[5]*M[0]*M[4]-kb[1]*c[5]*M[1]*M[0]*M[4]+
			kf[6]*rho*c[0]*c[4]*M[1]*M[5])/M[2]/M[1]/M[5]/M[4];
    dSwdU(NUM_VAR,NUM_VAR+3) += -kf[1]*rho*c[0]/M[3];
    dSwdU(NUM_VAR,NUM_VAR+4) += -1/M[0]*rho*(kb[3]*c[0]*M[0]*M[2]*M[5]*M[1]+kf[6]*rho*c[0]*
			       c[2]*M[0]*M[5]*M[1]-kb[6]*c[5]*M[0]*M[0]*M[2]*M[1]+2.0*kf[7]*rho*c[0]*c[0]*M[2]
			       *M[5]*M[1]-2.0*kb[7]*c[1]*M[0]*M[0]*M[2]*M[5])/M[4]/M[2]/M[5]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+5) += rho*(-kb[0]*c[0]*M[5]*M[2]*M[4]+kb[1]*c[2]*M[5]*M[0]*M[4]+
		       2.0*kf[3]*c[5]*M[0]*M[2]*M[4]+kb[6]*c[4]*M[5]*M[0]*M[2])/(M[5]*M[5])/M[2]/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR) += M[1]*rho*(kb[0]*c[5]*M[0]*M[4]+2.0*kf[7]*rho*c[0]*c[4]*M[5])/M[5]/(M[0]*M[0])/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -rho*(kf[0]*c[2]*M[4]+kb[7]*c[4]*M[2])/M[2]/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+2) += -kf[0]*rho/M[2]*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+4) += rho*(kf[7]*rho*c[0]*c[0]*M[1]-kb[7]*c[1]*M[0]*M[0])/(M[0]*M[0])/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+5) += M[1]*kb[0]*rho/M[5]*c[0]/M[0];
   
    dSwdU(NUM_VAR+2,NUM_VAR) += -rho*(-kb[0]*c[5]*M[3]*M[4]*M[2]-kf[1]*c[3]*M[5]*M[4]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[3])/M[0]/M[5]/M[3]/M[4];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += -kf[0]*rho*c[2]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -1/M[2]*rho*(kf[0]*c[1]*M[2]*M[5]*M[4]*M[0]+kb[1]*c[5]*M[2]
			       *M[1]*M[4]*M[0]+kb[2]*c[4]*M[2]*M[1]*M[5]*M[0]+4.0*kf[4]*rho*c[2]*c[4]*M[1]*M
			       [5]*M[0]+kf[5]*rho*c[5]*c[4]*M[2]*M[1]*M[0]+kf[6]*rho*c[0]*c[4]*M[2]*M[1]*M[5])/M[1]/M[5]/M[4]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += M[2]*rho*(kf[1]*c[0]*M[5]*M[4]+kf[2]*c[5]*M[0]*M[4]+2.0*kb
			    [4]*c[4]*M[0]*M[5])/M[0]/M[3]/M[5]/M[4];
    dSwdU(NUM_VAR+2,NUM_VAR+4) += -1/M[2]*rho*(kb[2]*c[2]*M[2]*M[4]*M[3]*M[5]*M[0]+2.0*kf[4]*
			       rho*c[2]*c[2]*M[4]*M[3]*M[5]*M[0]-2.0*kb[4]*c[3]*M[2]*M[2]*M[4]*M[5]*M[0]+kf[5]
			       *rho*c[2]*c[5]*M[2]*M[4]*M[3]*M[0]-2.0*kb[5]*c[4]*M[2]*M[2]*M[3]*M[5]*M[0]+kf
			       [6]*rho*c[0]*c[2]*M[2]*M[4]*M[3]*M[5]-kb[6]*c[5]*M[2]*M[2]*M[4]*M[3]*M[0])/(M[4]*M[4])/M[3]/M[5]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+5) += -rho*(-kb[0]*c[0]*M[3]*M[4]*M[2]+kb[1]*c[2]*M[0]*M[3]*M[4]-
			kf[2]*c[3]*M[0]*M[4]*M[2]+kf[5]*rho*c[2]*c[4]*M[0]*M[3]-kb[6]*c[4]*M[0]*M[3]*M
			  [2])/M[0]/M[5]/M[3]/M[4]; 
    dSwdU(NUM_VAR+3,NUM_VAR) += -kf[1]*rho/M[0]*c[3];
    dSwdU(NUM_VAR+3,NUM_VAR+2) += M[3]*rho*(kb[1]*c[5]*M[2]*M[4]+kb[2]*c[4]*M[5]*M[2]+2.0*kf
			    [4]*rho*c[2]*c[4]*M[5])/M[5]/(M[2]*M[2])/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+3) += -rho*(kf[1]*c[0]*M[5]*M[4]+kf[2]*c[5]*M[0]*M[4]+kb[4]*c[4]*
			M[0]*M[5])/M[0]/M[5]/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+4) += rho*(kb[2]*c[2]*M[2]*M[3]+kf[4]*rho*c[2]*c[2]*M[3]-kb[4]*c
		       [3]*M[2]*M[2])/(M[2]*M[2])/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+5) += -rho*(-kb[1]*c[2]*M[3]+kf[2]*c[3]*M[2])/M[5]/M[2];
   
    dSwdU(NUM_VAR+4,NUM_VAR) += -kb[3]*rho*c[4]/M[0];
    dSwdU(NUM_VAR+4,NUM_VAR+2) += rho*c[4]*(-kb[2]*M[5]+kf[5]*rho*c[5])/M[2]/M[5];
    dSwdU(NUM_VAR+4,NUM_VAR+3) += M[4]*kf[2]*rho*c[5]/M[5]/M[3];
    dSwdU(NUM_VAR+4,NUM_VAR+4) += 1/M[4]*rho*(-kb[2]*c[2]*M[4]*M[0]*M[5]-kb[3]*c[0]*M[2]*M[4]
			      *M[5]+kf[5]*rho*c[2]*c[5]*M[4]*M[0]-2.0*kb[5]*c[4]*M[2]*M[0]*M[5])/M[2]/M[0]/M[5];
    dSwdU(NUM_VAR+4,NUM_VAR+5) += rho*(kf[2]*c[3]*M[5]*M[2]*M[4]+2.0*kf[3]*c[5]*M[3]*M[2]*M
		       [4]+kf[5]*rho*c[2]*c[4]*M[5]*M[3])/(M[5]*M[5])/M[3]/M[2];
    dSwdU(NUM_VAR+5,NUM_VAR) += rho*(-kb[0]*c[5]*M[3]*M[4]*M[2]+kf[1]*c[3]*M[5]*M[4]*M[2]+
		       2.0*kb[3]*c[4]*M[5]*M[3]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[3])/M[0]/M[3]/M[4]/M[2];
    dSwdU(NUM_VAR+5,NUM_VAR+1) += M[5]*kf[0]*rho*c[2]/M[2]/M[1];
    dSwdU(NUM_VAR+5,NUM_VAR+2) += -rho*(-kf[0]*c[1]*M[5]*M[0]*M[4]+kb[1]*c[5]*M[1]*M[0]*M[4]-
			kb[2]*c[4]*M[1]*M[5]*M[0]+kf[5]*rho*c[5]*c[4]*M[1]*M[0]-kf[6]*rho*c[0]*c[4]*M
			[1]*M[5])/M[2]/M[1]/M[0]/M[4];
    dSwdU(NUM_VAR+5,NUM_VAR+3) += -rho*(-kf[1]*c[0]*M[5]+kf[2]*c[5]*M[0])/M[0]/M[3];
    dSwdU(NUM_VAR+5,NUM_VAR+4) += -rho*(-kb[2]*c[2]*M[4]*M[0]*M[5]-2.0*kb[3]*c[0]*M[2]*M[4]*M
			[5]+kf[5]*rho*c[2]*c[5]*M[4]*M[0]-2.0*kb[5]*c[4]*M[2]*M[0]*M[5]-kf[6]*rho*c[0]*
			c[2]*M[4]*M[5]+kb[6]*c[5]*M[2]*M[4]*M[0])/M[2]/(M[4]*M[4])/M[0];
    dSwdU(NUM_VAR+5,NUM_VAR+5) += -1/M[5]*rho*(kb[0]*c[0]*M[5]*M[2]*M[3]*M[4]+kb[1]*c[2]*M[5]
			       *M[0]*M[3]*M[4]+kf[2]*c[3]*M[5]*M[0]*M[2]*M[4]+4.0*kf[3]*c[5]*M[0]*M[2]*M[3]*M
			       [4]+kf[5]*rho*c[2]*c[4]*M[5]*M[0]*M[3]+kb[6]*c[4]*M[5]*M[0]*M[2]*M[3])/M[0]/M[2]/M[3]/M[4];    
      
    break;
    

    // 15 step CH4 based on GRI 2.11
  case CH4_15STEP_ARM2:
    Complex_Step_dSwdU(dSwdU, W);

//     //cout << "\n Complex Step dSwdU: " << dSwdU;
//     //dSwdU.zero();
//     //Finite_Difference_dSwdU(dSwdU, W, Flow_Type, Simple_Chemistry);
//     //cout << "\n Finite Difference dSwdU: " << dSwdU;
    break;


    // 15 step CH4 based on GRI 3
  case CH4_15STEP_ARM3:
    Complex_Step_dSwdU(dSwdU, W);
    break;


    //---------------------------------//
    //----- User Specified ------------//
    //---------------------------------//
  case USER:
    cerr<<"\nUser specified not set up yet";
    exit(1);
    break;
  default:
    //Do nothing (i.e. Jacobian = ZERO)
    break;
  };

    
} //end dSwdU





#endif // _LESPREMIXED2D_REACTIONS_INCLUDED

/**********************************************************************
 * HighTemp2DState.cc                                                 *
 *                                                                    *
 * Subroutines for 2D N-S HighTemp solution state classes.            *
 *                                                                    *
 **********************************************************************/

#include "HighTemp2DState.h"

/**********************************************************************
 * HighTemp2D_pState -- Create storage and assign gas constants.  *
 **********************************************************************/
int HighTemp2D_pState::eos_type = EOS_TGAS;
double HighTemp2D_pState::g    = GAMMA_AIR;
double HighTemp2D_pState::gm1  = GAMMA_AIR-ONE;
double HighTemp2D_pState::gm1i = ONE/(GAMMA_AIR-ONE);
//double HighTemp2D_pState::g    = 0;   
//double HighTemp2D_pState::gm1  = 0;
//double HighTemp2D_pState::gm1i = 0;
double HighTemp2D_pState::R    = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
double HighTemp2D_pState::cp   = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double HighTemp2D_pState::cv   = (R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double HighTemp2D_pState::v1   = AIR_c1;
double HighTemp2D_pState::v2   = AIR_c2;
double HighTemp2D_pState::v3   = AIR_c3;
double HighTemp2D_pState::v4   = AIR_c4;
double HighTemp2D_pState::v5   = AIR_c5;
double HighTemp2D_pState::Mref   = 0.5;
int HighTemp2D_pState::flow_type = FLOWTYPE_LAMINAR;

/**********************************************************************
 * HighTemp2D_cState -- Create storage and assign gas constants.  *
 **********************************************************************/
int HighTemp2D_cState::eos_type = EOS_TGAS;
double HighTemp2D_cState::g    = GAMMA_AIR;
double HighTemp2D_cState::gm1  = GAMMA_AIR-ONE;
double HighTemp2D_cState::gm1i = ONE/(GAMMA_AIR-ONE);
//double HighTemp2D_cState::g    = 0;
//double HighTemp2D_cState::gm1  = 0;
//double HighTemp2D_cState::gm1i = 0;
double HighTemp2D_cState::R    = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
double HighTemp2D_cState::cp   = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double HighTemp2D_cState::cv   = (R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double HighTemp2D_cState::v1   = AIR_c1;
double HighTemp2D_cState::v2   = AIR_c2;
double HighTemp2D_cState::v3   = AIR_c3;
double HighTemp2D_cState::v4   = AIR_c4;
double HighTemp2D_cState::v5   = AIR_c5;
int HighTemp2D_cState::flow_type = FLOWTYPE_LAMINAR;

/**********************************************************************
 * HighTemp2D_pState -- Create storage and assign turbulent       *
 *                          static variables.                         *
 **********************************************************************/
// Turbulent boundary-layer constants:
double HighTemp2D_pState::yplus_o = 10.0;
double HighTemp2D_pState::C = 5.0;
double HighTemp2D_pState::von_karman = 0.41;
double HighTemp2D_pState::yplus_sublayer = 5.0;
double HighTemp2D_pState::yplus_buffer_layer = 30.0;
double HighTemp2D_pState::yplus_outer_layer = 100.0;
// k-omega closure coefficients:
double HighTemp2D_pState::PrT = 0.90;
double HighTemp2D_pState::Cmu = 0.090;
double HighTemp2D_pState::beta_k_o = 0.09;
double HighTemp2D_pState::beta_omega_o = 0.072;
double HighTemp2D_pState::sigma_k = 0.50;
double HighTemp2D_pState::sigma_omega = 0.50;
double HighTemp2D_pState::alpha = 0.52;
double HighTemp2D_pState::xi = 1.50;
double HighTemp2D_pState::Mto = 0.25;

/**********************************************************************
 * HighTemp2D_cState -- Create storage and assign turbulent       *
 *                          static variables.                         *
 **********************************************************************/
// Turbulent boundary-layer constants:
double HighTemp2D_cState::yplus_o = 10.0;
double HighTemp2D_cState::C = 5.0;
double HighTemp2D_cState::von_karman = 0.41;
double HighTemp2D_cState::yplus_sublayer = 5.0;
double HighTemp2D_cState::yplus_buffer_layer = 30.0;
double HighTemp2D_cState::yplus_outer_layer = 100.0;
// k-omega coefficients:
double HighTemp2D_cState::PrT = 0.90;
double HighTemp2D_cState::Cmu = 0.090;
double HighTemp2D_cState::beta_k_o = 0.09;
double HighTemp2D_cState::beta_omega_o = 0.072;
double HighTemp2D_cState::sigma_k = 0.50;
double HighTemp2D_cState::sigma_omega = 0.50;
double HighTemp2D_cState::alpha = 0.52;
double HighTemp2D_cState::xi = 1.50;
double HighTemp2D_cState::Mto = 0.25;

/**********************************************************************
 * HighTemp2D_pState -- Create storage and assign propellant      *
 *                          constants.                                *
 **********************************************************************/
double HighTemp2D_pState::rhos = RHOS_APHTPB;
double HighTemp2D_pState::n    = N_APHTPB;
double HighTemp2D_pState::beta = BETA_APHTPB;
double HighTemp2D_pState::Tf   = TF_APHTPB;
double HighTemp2D_pState::Ts   = TS_APHTPB;

/**********************************************************************
 * HighTemp2D_cState -- Create storage and assign propellant      *
 *                          constants.                                *
 **********************************************************************/
double HighTemp2D_cState::rhos = RHOS_APHTPB;
double HighTemp2D_cState::n    = N_APHTPB;
double HighTemp2D_cState::beta = BETA_APHTPB;
double HighTemp2D_cState::Tf   = TF_APHTPB;
double HighTemp2D_cState::Ts   = TS_APHTPB;

/**********************************************************************
 * Routine: Riemann (Exact Riemann solver, x-direction)               *
 *                                                                    *
 * This function uses a Newton-Raphson interative procedure to obtain *
 * the exact solution to the Riemann problem for the 2D N-S HighTemp  *
 * equations in the x-direction, returning the intermediate state     *
 * variables along the ray x/t=0.  See Gottlieb and Groth (1987).     *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState Riemann(const HighTemp2D_pState &Wl,
			      const HighTemp2D_pState &Wr) {

  int number_of_iterations;
  
  double al, ar, CL, CR, Z;
  double dml, dmr, vm, vml, vmr, pm, aml, amr;
  double msl, pml, dpmldum, msr, pmr, dpmrdum;
  double vsl, vhl, vtl, vsr, vhr, vtr;

  HighTemp2D_pState W; W.Vacuum();

  // Determine the left and right state sound speeds.
  al = Wl.a(); ar = Wr.a();

  // Compute the left and right state Riemann invariants.
  CL = Wl.v.x + TWO*al/Wl.gm1;
  CR = Wr.v.x - TWO*ar/Wr.gm1;

  // Check for vacuum state.
  if (CL-CR <= ZERO) return W;

  // Make an initial estimate of the intermediate state flow velocity to
  // begin the Newton-Raphson iterative solution procedure.  The initial
  // guess tate velocity is made based on isentropic flow theory.
  Z = (ar/al)*pow((Wl.p/Wr.p),HALF*Wl.gm1/Wl.g);
  vm = (CL*Z + CR)/(ONE + Z);
  
  // In the case that two rarefaction waves are present, then an exact
  // solution has been found and the iterative procedure is not 
  // required.  Check for this.
  if (vm >= Wl.v.x && vm <= Wr.v.x) {
    if (vm >= ZERO) {
      aml = al - HALF*Wl.gm1*(vm - Wl.v.x);
      pm = Wl.p*pow((aml/al),TWO*Wl.g/Wl.gm1);
      vhl = Wl.v.x - al;
      vtl = vm - aml;
      if (vhl >= ZERO) {
	return Wl;
      } else if (vtl <= ZERO) {
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(aml);
	W = HighTemp2D_pState(dml,vm,vml,pm);
	return W;
      } else {
	vm = (Wl.gm1*Wl.v.x + TWO*al)/(Wl.g + ONE);
	pm = Wl.p*pow((vm/al),TWO*Wl.g/Wl.gm1);
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(vm);
	W = HighTemp2D_pState(dml,vm,vml,pm);
	return W;
      }
    } else {
      amr = ar + HALF*Wr.gm1*(vm - Wr.v.x);
      pm = Wr.p*pow((amr/ar),TWO*Wr.g/Wr.gm1);
      vhr = Wr.v.x+ar;
      vtr = vm + amr;
      if (vhr <= ZERO) {
	return Wr;
      } else if (vtr >= ZERO) {
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(amr);
	W = HighTemp2D_pState(dmr,vm,vmr,pm);
 	return W;
      } else {
	vm = (Wr.gm1*Wr.v.x - TWO*ar)/(Wr.g + ONE);
	pm = Wr.p*pow((-vm/ar),TWO*Wr.g/Wr.gm1);
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(vm);
	W = HighTemp2D_pState(dmr,vm,vmr,pm);
	return W;
      }
    }
  }
  
  // Perform the Newton-Raphson iterative procedure and solve for
  // the velocity in the intermediate state.  During this iterative
  // process the pressure in the intermediate state is also found.
  number_of_iterations = 0;
  while (1) {
    // Update the iteration counter.
    number_of_iterations = number_of_iterations + 1;
    
    // Determine solution changes for left wave.
    if (vm < Wl.v.x) {
      msl = (Wl.g+ONE)*(vm - Wl.v.x)/(FOUR*al);
      msl = msl-sqrt(ONE + sqr(msl));
      pml = Wl.p*(ONE + Wl.g*(vm - Wl.v.x)*msl/al);
      dpmldum = TWO*Wl.g*Wl.p*cube(msl)/(al*(ONE + sqr(msl)));
    } else {
      aml = al - HALF*Wl.gm1*(vm - Wl.v.x);
      pml = Wl.p*pow((aml/al),TWO*Wl.g/Wl.gm1);
      dpmldum = -Wl.g*pml/aml;
    }
	  
    // Determine solution changes for right wave.
    if (vm > Wr.v.x) {
      msr = (Wr.g+ONE)*(vm - Wr.v.x)/(FOUR*ar);
      msr = msr+sqrt(ONE + sqr(msr));
      pmr = Wr.p*(ONE+Wr.g*(vm - Wr.v.x)*msr/ar);
      dpmrdum = TWO*Wr.g*Wr.p*cube(msr)/(ar*(ONE+sqr(msr)));
    } else {
      amr = ar+HALF*Wr.gm1*(vm - Wr.v.x);
      pmr = Wr.p*pow((amr/ar),TWO*Wr.g/Wr.gm1);
      dpmrdum = Wr.g*pmr/amr;
    }
    
    // Check for convergence (i.e., pml=pmr).
    if (fabs(ONE-pml/pmr) <= TOLER) break;
    
    // Compute next estimate for the intermediate state velocity, vm.
    vm = vm-(pml-pmr)/(dpmldum-dpmrdum);
    
  }
  
  pm = HALF*(pml+pmr);
  
  // Return the intermediate state solution.
  /*
  if (vm >= ZERO) {
    if (vm < Wl.v.x) {
      aml = al*sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p)/
	            ((Wl.g+ONE)+Wl.gm1*Wl.p/pm));
      vsl = Wl.v.x+msl*al;
      if (vsl >= ZERO) {
	return HighTemp2D_pState(Wl);                
      } else {
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(aml);
	return HighTemp2D_pState(dml,vm,vml,pm);
      }
    } else {
      vhl = Wl.v.x-al;
      vtl = vm-aml;
      if (vhl >= ZERO) {
	return Wl;
      } else if (vtl <= ZERO) {
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(aml);
	W = HighTemp2D_pState(dml,vm,vml,pm);
	return W;
      } else {
	vm = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
	pm = Wl.p*pow((vm/al),TWO*Wl.g/Wl.gm1);
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(vm);
	W = HighTemp2D_pState(dml,vm,vml,pm);
	return W;
      }
    }
  } else {
  */
    if (vm > Wr.v.x) {
      amr = ar*sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p)/
		    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm));
      vsr = Wr.v.x+msr*ar;
      if (vsr <= ZERO) {
	return Wr;
      } else {
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(amr);
	W = HighTemp2D_pState(dmr,vm,vmr,pm);
	return W;
      }
    } else {
      vhr = Wr.v.x + ar;
      vtr = vm + amr;
      if (vhr <= ZERO) {
	return HighTemp2D_pState(Wr);
      } else if (vtr >= ZERO) {
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(amr);
	W = HighTemp2D_pState(dmr,vm,vmr,pm);
	return W;
      } else {
	vm = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
	pm = Wr.p*pow((-vm/ar),TWO*Wr.g/Wr.gm1);
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(vm);
	W = HighTemp2D_pState(dmr,vm,vmr,pm);
	return W;
      }
    }
    //  }

}

/**********************************************************************
 * Routine: RoeAverage (Roe Averages)                                 *
 *                                                                    *
 * This function returns the Roe-averaged primitive solution state    *
 * given left and right primitive solution variables. See Roe (1981). *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState RoeAverage(const HighTemp2D_pState &Wl, const HighTemp2D_pState &Wr) {
  
  double hl, hr, srhol, srhor, aa2, ha;
  HighTemp2D_pState Wa; Wa.rho = ZERO; Wa.p = ZERO;//Wa.Vacuum();

  // Determine the left and right state specific enthalpies and square
  // roots of the density.
  hl    = Wl.h();
  hr    = Wr.h();
  srhol = sqrt(Wl.rho);
  srhor = sqrt(Wr.rho);
  
  // Determine the appropriate Roe averages.
  Wa.rho   = srhol*srhor;
  Wa.v.x   = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
  Wa.v.y   = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
  Wa.k     = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
  Wa.omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
  ha  = (srhol*hl + srhor*hr)/(srhol+srhor);

  //added b/c gamma not defined anymore :) - just for testing
  // Wa.gm1 = GAMMA_AIR-ONE;  
  // Wa.g = GAMMA_AIR;

  aa2 = Wa.gm1*(ha - HALF*(sqr(Wa.v.x) + sqr(Wa.v.y)) - Wa.k);
  Wa.p = Wa.rho*aa2/Wa.g;

  double el, er, ea;
  el = Wl.e();
  er = Wr.e();
  ea = (srhol*el + srhor*er)/(srhol+srhor);
  //cout<<"RoeAverage Wa.rho = "<<Wa.rho<<" Wa.p = "<<Wa.p<<" aa2 = "<<aa2<<endl;

  // Return the Roe-averged state.
  return Wa;

}

/**********************************************************************
 * Routine: Translate                                                 *
 *                                                                    *
 * This function returns the solution in a stationary frame.          *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState Translate(const HighTemp2D_pState &W, const Vector2D &V) {

  return HighTemp2D_pState(W.rho,W.v.x-V.x,W.v.y-V.y,W.p,W.k,W.omega);

}

/**********************************************************************
 * Routine: Reflect                                                   *
 *                                                                    *
 * This function returns the reflected solution state in a given      *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState Reflect(const HighTemp2D_pState &W,
			      const Vector2D &norm_dir) {

  HighTemp2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Reflect the gas-phase normal velocity in the rotated frame.
  Wr.v.x = - Wr.v.x;

  // Rotate back to the original Cartesian reference frame.
  Wn.Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the reflected state.
  return Wn;

}

/**********************************************************************
 * Routine: Mirror                                                    *
 *                                                                    *
 * This function returns the mirrored solution state in a given       *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState Mirror(const HighTemp2D_pState &W, const Vector2D &norm_dir) {
  
  HighTemp2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Mirror the normal velocity in the rotated frame.
  Wr.v.x = -Wr.v.x;
  Wr.v.y = -Wr.v.y;

  // Rotate back to the original Cartesian reference frame.
  Wn.Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the mirrored state.
  return Wn;
  
}

/**********************************************************************
 * Routine: WallViscousHeatFlux                                       *
 *                                                                    *
 * This function returns the adiabatic no-slip solution state in a    *
 * given direction given the primitive solution variables and the     *
 * unit normal vector in the direction of interest.                   *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState WallViscousHeatFlux(const HighTemp2D_pState &W,
					  const Vector2D &norm_dir) {

  HighTemp2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Use the Roe-averaged value to find the correct tangential velocity.
  Wr.v.x = - Wr.v.x;
  Wr.v.y = - Wr.v.y;

  // Apply the adiabatic viscous wall boundary condition to the
  // turbulent variables.
  Wr.k = ZERO;

  // Rotate back to the original Cartesian reference frame.
  Wn.Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the moving wall state.
  return Wn;
//   return MovingWallHeatFlux(W,norm_dir,ZERO);

}

/**********************************************************************
 * Routine: WallViscousIsothermal                                     *
 *                                                                    *
 * This function returns the no-slip solution state with a fixed wall *
 * temperature in a given direction given the primitive solution      *
 * variables and the unit normal vector in the direction of interest. *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState WallViscousIsothermal(const HighTemp2D_pState &W,
				     const Vector2D &norm_dir,
				     const double &Twall) {

  return MovingWallIsothermal(W,norm_dir,ZERO,Twall);

}

/**********************************************************************
 * Routine: MovingWallHeatFlux                                        *
 *                                                                    *
 * This function returns the adiabatic moving wall solution state in  *
 * a given direction given the primitive solution variables and the   *
 * unit normal vector in the direction of interest.                   *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState MovingWallHeatFlux(const HighTemp2D_pState &W,
					 const Vector2D &norm_dir,
					 const double &Vwall) {

  HighTemp2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Use the Roe-averaged value to find the correct tangential velocity.
  Wr.v.x = - Wr.v.x;
  Wr.v.y = - TWO*Vwall - Wr.v.y;

  // Apply the adiabatic viscous wall boundary condition to the
  // turbulent variables.
  Wr.k = ZERO;

  // Rotate back to the original Cartesian reference frame.
  Wn.Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the moving wall state.
  return Wn;

}

/**********************************************************************
 * Routine: MovingWallIsothermal                                      *
 *                                                                    *
 * This function returns the moving wall solution state with a fixed  *
 * wall temperature in a given direction given the primitive solution *
 * variables and the  unit normal vector in the direction of          *
 * interest.                                                          *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState MovingWallIsothermal(const HighTemp2D_pState &W,
					   const Vector2D &norm_dir,
					   const double &Vwall,
					   const double &Twall) {

  HighTemp2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Use the Roe-averaged value to find the correct tangential velocity.
  Wr.v.x = - Wr.v.x;
  Wr.v.y = - TWO*Vwall - Wr.v.y;

  // Apply the specified wall temperature(modified for HighTemp2D cases)
  switch(W.eos_type){
    case EOS_TGAS:
      Wr.rho = Tgas_rho(W.p,Twall);
      break;
    case EOS_IDEAL:
      Wr.rho = W.p/(W.R*Twall);
      break;
  }

  // Apply the isothermal viscous wall boundary condition to the
  // turbulent variables.
  Wr.k = ZERO;

  // Rotate back to the original Cartesian reference frame.
  Wn.Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the moving wall state.
  return Wn;
       
}

/**********************************************************************
 * Routine: BurningSurface                                            *
 *                                                                    *
 * This function returns the star state solution in a given direction *
 * for a burning surface given the primitive solution variables and   *
 * the unit normal vector in the direction of interest.  Burning rate *
 * based on Saint Robert's pressure dependent burning law.            *
 *                                                                    *
 **********************************************************************/

HighTemp2D_pState BurningSurface(const HighTemp2D_pState &W,
				     const Vector2D &norm_dir) {

  double g = W.g, R = W.R, vxr, vyr, pr, ar; 
  double rhobs, vxbs, vxbsp, rbs = ZERO, rbsp, vybs, pbs, pbso, ars;
  double cos_angle, sin_angle, vx, vy, C1, C2;
  HighTemp2D_pState Wbs;
  int iteration = 0;

  // Determine the direction cosine's for the frame rotation.
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  vxr =   W.v.x*cos_angle + W.v.y*sin_angle;
  vyr = - W.v.x*sin_angle + W.v.y*cos_angle;
  pr  = W.p;
  ar  = W.a();
  vxr = - vxr; // Required to correct for burning direction.

  // Compute the star state at the burning surface.
  // initial guess.
  pbs = pr;
  iteration = 0;
  do {
    pbso = pbs;
    if (pr > pbso) {
      // Rarefaction wave.
      ars   = ar*pow(pbso/pr,(g-ONE)/(TWO*g));
      vxbs  = vxr + (TWO/(g-ONE))*(ars-ar);
      vxbsp = ars/(g*pbso);
    } else if (pr <= pbso) {
      // Shock wave.
      C1    = ((g+ONE)/(TWO*g))*pbso/pr + ((g-ONE)/(TWO*g));
      C2    = ONE - ((g+ONE)/(FOUR*g))*(pbso/pr - ONE)/C1;
      vxbs  = vxr + ar*(pbso/pr - ONE)/(g*sqrt(C1));
      vxbsp = (ar*C2/(g*pr))/sqrt(C1);
    }
    rbs   = W.beta*pow(pbso,W.n);
    rbsp  = W.n*rbs/pbso;
    rhobs = pbso/(R*W.Tf);
    pbs   = pbso - (W.rhos*rbs - rhobs*vxbs)/(W.rhos*rbsp - rhobs*(vxbs/pbso + vxbsp));
    if (pbs <= ZERO) pbs = HALF*pbso;
    iteration++;
  } while (fabs(ONE - (pbso/pbs)) > TOLER/HUNDRED);
  vybs = ZERO;
  vxbs = -vxbs; // Required to correct for burning direction.

  // Rotate back to the original reference frame.
  Wbs.rho   = rhobs;
  Wbs.v.x   = vxbs*cos_angle - vybs*sin_angle;
  Wbs.v.y   = vxbs*sin_angle + vybs*cos_angle;
  Wbs.p     = pbs;
  if (Wbs.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    double sigmav = 0.035;
    double lw = 0.000200;
    Wbs.k     = sqr(sigmav)*sqr(vxbs);
    Wbs.omega = pow(W.Cmu,-0.75)/(lw*Wbs.beta_k_o);
  }

  // Return the burning surface state.
  return Wbs;

}

/**********************************************************************
 * Routine: RinglebFlow                                               *
 *                                                                    *
 * This function returns the exact solution to Ringleb's flow for the *
 * location X.                                                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState RinglebFlow(const HighTemp2D_pState &Wdum,
				  const Vector2D &X) {
  double q, k;
  return RinglebFlow(Wdum,X,q,k);
}

HighTemp2D_pState RinglebFlow(const HighTemp2D_pState &Wdum,
				  const Vector2D &X,
				  double &q, double &k) {

  HighTemp2D_pState W; W.Vacuum();
  double sin_theta, cos_theta, theta;
  double f_a, f_ab;
  double J, J_a, J_ab;
  double rho_a, rho_ab;
  double q_a, q_ab;
  double c, c_a = 0.70, c_b = 1.00, c_ab;
  double g = GAMMA_AIR, po = PRESSURE_STDATM, rhoo = DENSITY_STDATM;

  // Use bisection method to solve for the sound speed, c.
  while (fabs(c_a - c_b) > TOLER*TOLER) {
    // Determine f_a.
    rho_a = pow(c_a,TWO/(g-ONE));
    J_a = ONE/c_a + ONE/(THREE*c_a*c_a*c_a) + ONE/(FIVE*c_a*c_a*c_a*c_a*c_a) - HALF*log((ONE+c_a)/(ONE-c_a));
    q_a = sqrt((TWO/(g-ONE))*(ONE-c_a*c_a));
    f_a = (X.x + HALF*J_a)*(X.x + HALF*J_a) + X.y*X.y - ONE/(FOUR*rho_a*rho_a*q_a*q_a*q_a*q_a);
    // Determine f_ab.
    c_ab = HALF*(c_a + c_b);
    rho_ab = pow(c_ab,TWO/(g-ONE));
    J_ab = ONE/c_ab + ONE/(THREE*c_ab*c_ab*c_ab) + ONE/(FIVE*c_ab*c_ab*c_ab*c_ab*c_ab) - HALF*log((ONE+c_ab)/(ONE-c_ab));
    q_ab = sqrt((TWO/(g-ONE))*(ONE-c_ab*c_ab));
    f_ab = (X.x + HALF*J_ab)*(X.x + HALF*J_ab) + X.y*X.y - ONE/(FOUR*rho_ab*rho_ab*q_ab*q_ab*q_ab*q_ab);
    // Update the bounds of the bisection search as required.
    if (f_a*f_ab <= ZERO) c_b = HALF*(c_a + c_b);
    else c_a = HALF*(c_a + c_b);
  }

  // Final sound speed, density, and total velocity (speed).
  c = HALF*(c_a + c_b);
  q = sqrt((TWO/(g-ONE))*(ONE-c*c));
  W.rho = pow(c,TWO/(g-ONE));
  J = ONE/c + ONE/(THREE*c*c*c) + ONE/(FIVE*c*c*c*c*c) - HALF*log((ONE+c)/(ONE-c));
  k = sqrt(TWO/(TWO*W.rho*(X.x+HALF*J)+ONE/(q*q)));
  //if (k > 5.0/3.0) cout << "k = " << k << " > 5/3 @ " << X << endl;
  sin_theta = max(ZERO,min(ONE,q/k));
  theta = TWO*PI-asin(sin_theta);
  sin_theta = sin(theta);
  cos_theta = cos(theta);
  W.rho = rhoo*W.rho;
  W.v.x = sqrt(g*po/rhoo)*q*cos_theta;
  if (X.y < ZERO) W.v.x = -ONE*W.v.x;
  W.v.y = sqrt(g*po/rhoo)*q*sin_theta;
  W.p   = po*(W.rho/rhoo)*c*c;

  // Return W state.
  return W;

}

/**********************************************************************
 * Routine: RinglebFlow                                               *
 *                                                                    *
 * This function returns the average exact solution to Ringleb's flow *
 * for a quadrilateral cell defined by points X1, X2, X3, and X4      *
 * which are defined in a counter-clockwise direction starting in the *
 * south-west corner.  The average solution is determined using a     *
 * ?????th-order numerical quadrature.                                *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState RinglebFlowAverageState(const HighTemp2D_pState &Wdum,
					      const Vector2D &Y1,
					      const Vector2D &Y2,
					      const Vector2D &Y3,
					      const Vector2D &Y4) {

  Vector2D X, X1, X2, X3, X4;
  double epsilon1 = -ONE, epsilon2 = ONE, epsilon3 = -ONE, epsilon4 = ONE;
  double eta1 = -ONE, eta2 = -ONE, eta3 =  ONE, eta4 = ONE;
  double N1, N2, N3, N4;
  double w[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
  double epsilon[3] = {-sqrt(15.0)/5.0, ZERO, sqrt(15.0)/5.0};
  double eta[3] = {-sqrt(15.0)/5.0, ZERO, sqrt(15.0)/5.0};
  double u[7] = {0.225000000000000,
		 0.125939180544827, 0.125939180544827,
		 0.125939180544827, 0.132394152788506,
		 0.132394152788506, 0.132394152788506};
  double upsilon1[7] = {1.0/3.0,
			0.797426985353087, 0.101286507323456,
			0.101286507323456, 0.059715871789770,
			0.470142064105115, 0.470142064105115};
  double upsilon2[7] = {1.0/3.0,
			0.101286507323456, 0.797426985353087,
			0.101286507323456, 0.470142064105115,
			0.059715871789770, 0.470142064105115};
  double upsilon3[7] = {1.0/3.0,
			0.101286507323456, 0.101286507323456,
			0.797426985353087, 0.470142064105115,
			0.470142064105115, 0.059715871789770};
  HighTemp2D_pState W; W.Vacuum();

  if (abs(Y1-Y2) > NANO && abs(Y1-Y4) > NANO &&
      abs(Y2-Y3) > NANO && abs(Y3-Y4) > NANO) {

    // Perform the quintic numerical Gauss-quadrature over the quadrilateral.
    X1 = Y1; X2 = Y2; X3 = Y3; X4 = Y4;

    // Get the value of 'f' at each point.
    for (int jj = 0; jj < 3; jj++) {
      for (int ii = 0; ii < 3; ii++) {
	// Set basis functions.
	N1 = 0.25*(ONE + epsilon[ii]*epsilon1)*(ONE + eta[jj]*eta1);
	N2 = 0.25*(ONE + epsilon[ii]*epsilon2)*(ONE + eta[jj]*eta2);
	N3 = 0.25*(ONE + epsilon[ii]*epsilon3)*(ONE + eta[jj]*eta3);
	N4 = 0.25*(ONE + epsilon[ii]*epsilon4)*(ONE + eta[jj]*eta4);
	// Get point X.
	X = X1*N1 + X2*N2 + X3*N3 + X4*N4;
	// Determine the value of 'f'.
	W += w[ii]*w[jj]*RinglebFlow(W,X)/FOUR;
      }
    }

  } else {

    // Perform the quintic numerical Gauss-quadrature over the triangle.
    if (abs(Y1-Y2) < NANO) { X1 = Y1; X2 = Y3; X3 = Y4; }
    else if (abs(Y1-Y4) < NANO) { X1 = Y1; X2 = Y2; X3 = Y3; }
    else if (abs(Y2-Y3) < NANO) { X1 = Y1; X2 = Y2; X3 = Y4; }
    else if (abs(Y3-Y4) < NANO) { X1 = Y1; X2 = Y2; X3 = Y3; }

    double A = 0.5*((X1^X2) + (X2^X3) + (X3^X1));
    double d1 = X2.x*X3.y - X3.x*X2.y;
    double d2 = X3.x*X1.y - X1.x*X3.y;
    double d3 = X1.x*X2.y - X2.x*X1.y;
    double a1 = X3.x - X2.x;
    double a2 = X1.x - X3.x;
    double a3 = X2.x - X1.x;
    double b1 = X2.y - X3.y;
    double b2 = X3.y - X1.y;
    double b3 = X1.y - X2.y;
    double det = d1*b2*a3-d1*a2*b3-d2*b1*a3+d2*a1*b3+d3*b1*a2-d3*a1*b2;

    for (int nn = 0; nn < 7; nn++) {
      X.x = TWO*A*((a2*d3-d2*a3)*upsilon1[nn]+(d1*a3-a1*d3)*upsilon2[nn]+(a1*d2-d1*a2)*upsilon3[nn])/det;
      X.y = TWO*A*((d2*b3-b2*d3)*upsilon1[nn]+(b1*d3-d1*b3)*upsilon2[nn]+(d1*b2-b1*d2)*upsilon3[nn])/det;
      W += u[nn]*RinglebFlow(W,X);
    }

  }

  // Return the average W state of the specified cell.
  return W;

}

/**********************************************************************
 * Routine: ViscousChannelFlow                                        *
 *                                                                    *
 * This function will return the exact laminar channel flow solution  *
 * (Couette or Poiseuille flows) given an (x,y)-coordinate where x =  *
 * [0,length] and y = [0,height], an upper wall speed Vwall, and the  *
 * imposed pressure change, dp.  Use ViscousPipeFlow for axisymmetric *
 * problems.  The exact laminar solution can be used as an initial    *
 * solution for a turbulent channel flow.                             *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState ViscousChannelFlow(const HighTemp2D_pState &Wdum,
					 const Vector2D X,
					 const Vector2D Vwall,
					 const double dp,
					 const double length,
					 const double height) {
  HighTemp2D_pState W;
  // Compute the exact laminar channel solution.  Note that the
  // pressure, p, and density, rho, must be set before the x-direction
  // velocity component, v.x, since the v.x depends on the dynamic
  // viscosity, mu, which depends on the temperature, T.
  //  W.rho = DENSITY_STDATM;
  W.rho = Wdum.rho;
  W.p   = Wdum.p - dp*(ONE - X.x/length);
  //  W.p   = PRESSURE_STDATM - dp*(ONE - X.x/length);
  W.v.x = (HALF/W.mu())*dp*X.y*(X.y - height)/length + Vwall.x*X.y/height;
  W.v.y = ZERO;
  // Return W state.
  return W;
}

HighTemp2D_pState ViscousChannelFlowVelocity(const HighTemp2D_pState &Wdum,
						 const Vector2D X,
						 const Vector2D Vwall,
						 const double dp,
						 const double length,
						 const double height) {
  HighTemp2D_pState W;
  // Compute the exact laminar channel solution.  Note that the
  // pressure, p, and density, rho, must be set before the x-direction
  // velocity component, v.x, since the v.x depends on the dynamic
  // viscosity, mu, which depends on the temperature, T.
  W.rho = Wdum.rho;
  W.p   = Wdum.p;
  W.v.x = (HALF/W.mu())*dp*X.y*(X.y - height)/length + Vwall.x*X.y/height;
  W.v.y = ZERO;
  // Return W state.
  return W;
}

/**********************************************************************
 * Routine: ViscousPipeFlow                                           *
 *                                                                    *
 * This function will return the exact laminar flow solution for an   *
 * axisymmetric pipe flow, given the imposed pressure drop, dp, and   *
 * an (r,z)-coordinate where r = [0,radius] and z = [0,length].       *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState ViscousPipeFlow(const HighTemp2D_pState &Wdum,
				      const Vector2D X,
				      const double dp,
				      const double length,
				      const double radius) {

  HighTemp2D_pState W;
  // Compute the exact laminar pipe solution.  Note that the pressure,
  // p, and density, rho, must be set before the x-direction velocity 
  // component, v.x, since the v.x depends on the dynamic viscosity, 
  // mu, which depends on the temperature, T.
  W.rho = Wdum.rho;
  W.p   = Wdum.p - (dp*length)*(ONE - X.x/length);
  //W.rho = DENSITY_STDATM;
  // W.p   = PRESSURE_STDATM - (dp*length)*(ONE - X.x/length);
  W.v.x = (X.y*X.y - radius*radius)*dp*length/(FOUR*W.mu());
  W.v.y = ZERO;
  // Return W state.
  return W;

}

/**********************************************************************
 * Routine: TurbulentPipeFlow                                         *
 *                                                                    *
 * This function will return an approximate turbulent flow solution   *
 * for an axisymmetric pipe flow, given the imposed pressure drop,    *
 * dp, the dimensions of the pipe and an (r,z)-coordinate where r =   *
 * [0,radius] and z = [0,length].                                     *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState TurbulentPipeFlow(const HighTemp2D_pState &Wo,
					const Vector2D X,
					const double dp,
					const double length,
					const double radius) {

  HighTemp2D_pState W;
  double utau, tauw, y;
  // Compute the approximate turbulent pipe solution.  Note that the
  // pressure, p, and density, rho, must be set before the x-direction
  // velocity component, v.x, since the v.x depends on the dynamic
  // viscosity, mu, which depends on the temperature, T.
  y = max(ZERO,min(X.y,radius-TOLER));
  W.rho = Wo.rho;
  W.p   = Wo.p + X.x*dp/length;
  W.v.x = Wo.v.x*pow(1.0 - y/radius,0.133);
  W.v.y = ZERO;
  tauw = -HALF*radius*(dp/length);
  utau = sqrt(tauw/W.rho);
  W.k = sqr(utau)/sqrt(W.beta_k_o);
  W.omega = utau/(sqrt(W.beta_k_o)*W.von_karman*(radius-X.y));
  // Return W state.
  return W;

}

/**********************************************************************
 * Routine: FlatPlate                                                 *
 *                                                                    *
 * This function returns the exact solution for the flow over a flat  *
 * (adiabatic) plate (Blasius solution) at a given the position and   *
 * the freestream flow velocity.                                      *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState FlatPlate(const HighTemp2D_pState &Winf,
				const Vector2D &X,
				const double &plate_length,
				double &eta,
				double &f,
				double &fp,
				double &fpp) {

  HighTemp2D_pState W;
  double fo, n, dn, k1, k2, k3, k4;

  // Initialize variables.
  W.rho = Winf.rho; W.p = Winf.p;
  eta = ZERO; f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
  n = ZERO; dn = 0.00005;

  // Return upstream conditions before flat plate, including the leading edge.
  if (X.x < TOLER || X.y > HALF*plate_length) return Winf;
  //if (X.x < NANO || X.x > plate_length || X.y > plate_length) return Winf;

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(Winf.v.x/(X.x*Winf.nu()));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 100.0) {
    return Winf;
  } else if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    W.v.x = fp*Winf.v.x;
    W.v.y = sqrt(HALF*Winf.nu()*Winf.v.x/X.x)*(eta*fp-f);
    return W;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  while (n < eta) {

    // Store the solution at the start of the iteration.
    fo = f;

    // Increment f:
    k1 = dn*fp;
    k2 = dn*(fp + k1/2.0);
    k3 = dn*(fp + k2/2.0);
    k4 = dn*(fp + k3);
    f += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fp:
    k1 = dn*fpp;
    k2 = dn*(fpp + k1/2.0);
    k3 = dn*(fpp + k2/2.0);
    k4 = dn*(fpp + k3);
    fp += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fpp:
    k1 = -dn*fo*fpp/2.0;
    k2 = -dn*(fo + dn/2.0)*(fpp + k1/2.0)/2.0;
    k3 = -dn*(fo + dn/2.0)*(fpp + k2/2.0)/2.0;
    k4 = -dn*(fo + dn)*(fpp + k3)/2.0;
    fpp += k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Determine the increment dn:
    if (n + dn > eta) dn = eta - n;

    // Increment n:
    n += dn;

  }

  // Compute the velocity vector at point X.
  W.v.x = fp*Winf.v.x;
  W.v.y = sqrt(HALF*Winf.nu()*Winf.v.x/X.x)*(eta*fp-f);

  // Return the Blasius solution state.
  return W;

}
// HighTemp2D_pState FlatPlate(const HighTemp2D_pState &Winf,
// 				const Vector2D &X,
// 				const double &plate_length,
// 				double &eta,
// 				double &f,
// 				double &fp,
// 				double &fpp) {

//   HighTemp2D_pState W;
//   double fo, dn, k1, k2, k3, k4;

//   // Initialize variables.
//   W.rho = Winf.rho; W.p = Winf.p;
//   eta = ZERO; f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
//   dn = 0.00001;

//   // Return upstream conditions before flat plate, including the leading edge.
//   if (X.x < TOLER || X.y > HALF*plate_length) return Winf;
//   //if (X.x < NANO || X.x > plate_length || X.y > plate_length) return Winf;

//   // Determine the dimensionless similarity coordinate, eta:
//   eta = X.y*sqrt(Winf.v.x/(X.x*Winf.nu()));

//   // If eta is greater than 8.4, for the sake of expediency, use linear
//   // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
//   // given the tabulated value at 8.4 (note, the analytic solution is 
//   // linear in this region).
//   if (eta > 8.4) {
//     fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
//     W.v.x = fp*Winf.v.x;
//     W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/X.x)*(eta*fp-f);
//     return W;
//   }

//   // Compute the Blasius solution using a fourth-order Runge-Kutta method.
//   for (double n = ZERO; n < eta; n += dn) {

//     // Store the solution at the start of the iteration.
//     fo = f;

//     // Increment f:
//     k1 = dn*fp;
//     k2 = dn*(fp + k1/2.0);
//     k3 = dn*(fp + k2/2.0);
//     k4 = dn*(fp + k3);
//     f = f + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

//     // Increment fp:
//     k1 = dn*fpp;
//     k2 = dn*(fpp + k1/2.0);
//     k3 = dn*(fpp + k2/2.0);
//     k4 = dn*(fpp + k3);
//     fp = fp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

//     // Increment fpp:
//     k1 = -dn*fo*fpp/2.0;
//     k2 = -dn*(fo + dn/2.0)*(fpp + k1/2.0)/2.0;
//     k3 = -dn*(fo + dn/2.0)*(fpp + k2/2.0)/2.0;
//     k4 = -dn*(fo + dn)*(fpp + k3)/2.0;
//     fpp = fpp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

//   }

//   // Compute the velocity vector at point X.
//   W.v.x = fp*Winf.v.x;
//   W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/X.x)*(eta*fp-f);

//   // Return the Blasius solution state.
//   return W;

// }

/**********************************************************************
 * Routine: DrivenCavityFlow                                          *
 *                                                                    *
 * This function returns the initial state for the laminar flow in a  *
 * square cavity with a moving upper surface.  The initial state is   *
 * set so that the Mach number is 0.1 for cavity Reynolds numbers of  *
 * 100.0 or 400.0.  The density and pressure are set accordingly and  *
 * the gas is assumed to be air.                                      *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState DrivenCavityFlow(const HighTemp2D_pState &Wo,
				       const double &l,
				       const double &Re) {
  assert(fabs(Re-100.0) < TOLER || fabs(Re-400.0) < TOLER);
  assert(fabs(l-0.001) < TOLER);
  HighTemp2D_pState W;
  int n;
  double M, eta, c1, c2, c3, c4, c5, f, df, T, To, T1, T2, f1, f2;
  // Set the Mach number.
  M = 0.10;
  // Set the density and the predicted temperature.
  if (fabs(Re-100.0) < TOLER) {
    W.rho = 0.0265; T = 268.0; To = 270.0;
  } else if (fabs(Re-400.0) < TOLER) {
    W.rho = 0.1075; T = 335.0; To = 336.0;
  }
  // Determine the coefficients.
  c1 = AIR_c1; c2 = AIR_c2; c3 = AIR_c3; c4 = AIR_c4;
  c5 = TWO*W.rho*l*sqrt(W.g*W.R)*M/Re;
  // Determine the temperature.
  n = 0;
  while (fabs(T-To) > 1.0e-06) {
    To = T;
    f = c1*To - c5*(c2 + pow(To,c3) + c4/To);
    df = c1 - c5*(c3*pow(To,c3-ONE) - c4/(To*To));
    T = To - f/df;
    n++;
    if (n > 25) {
      cout << n << " " << T << " " << To << " " << f << " " << df << endl;
      break;
    }
  }
  // Determine the pressure.
  W.p = W.rho*W.R*T;
  // Set the velocity.
  W.v = Vector2D_ZERO;
  // Return the solution state.
  return W;
}

/**********************************************************************
 * Routine: BackwardFacingStep                                        *
 *                                                                    *
 * This function returns the initial conditions for the laminar flow  *
 * over a backward-facing step at a given position.  The initial      *
 * conditions are set so that the Mach number of the inflow (fully    *
 * developed pipe) is 0.2 for all three possible Reynolds numbers     *
 * (100.0, 389.0, 1000.0).  The density and pressure are set          *
 * accordingly.  The gas is assumed to be air.                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState BackwardFacingStep(const HighTemp2D_pState &Wo,
					 const Vector2D &X,
					 const double &h,
					 const double &ho,
					 const double &Re,
					 const double &M) {
  assert(fabs(Re-100.0) < TOLER || fabs(Re-389.0) < TOLER || fabs(Re-1000.0) < TOLER);
  assert(fabs(M-0.20) < TOLER);
  HighTemp2D_pState W;
  int n;
  double eta, c1, c2, c3, c4, c5, f, df, T, To, T1, T2, f1, f2;
  // Set the density.  Quite low eh.
  if (fabs(Re-100.0) < TOLER) W.rho = 0.0046;
  else if (fabs(Re-389.0) < TOLER) W.rho = 0.0104;
  else if (fabs(Re-1000.0) < TOLER) W.rho = 0.045;
  // Determine the temperature.
  c1 = AIR_c1; c2 = AIR_c2; c3 = AIR_c3; c4 = AIR_c4;
  c5 = TWO*W.rho*h*sqrt(W.g*W.R)*M/Re;
  T = 335.0; To = 336.0; n = 0;
  while (fabs(T-To) > 1.0e-06) {
    To = T;
    f = c1*To - c5*(c2 + pow(To,c3) + c4/To);
    df = c1 - c5*(c3*pow(To,c3-ONE) - c4/(To*To));
    T = To - f/df;
    n++;
    if (n > 25) {
      cout << n << " " << T << " " << To << " " << f << " " << df << endl;
      break;
    }
  }
  // Determine the pressure.
  W.p = W.rho*W.R*T;
  // Determine the x-component of the velocity.
  if (X.x > ZERO) eta = ZERO;
  else eta = ONE + TWO*(X.y - ho - h)/h;
  W.v.x = W.a()*M*(ONE - eta*eta);
  // Set the y-component of the velocity.
  W.v.y = ZERO;
  // Return the solution state.
  return W;
}

/**********************************************************************
 * Routine: ForwardFacingStep                                        *
 *                                                                    *
 * This function returns the initial conditions for the laminar flow  *
 * over a backward-facing step at a given position.  The initial      *
 * conditions are set so that the Mach number of the inflow (fully    *
 * developed pipe) is 0.2 for all three possible Reynolds numbers     *
 * (100.0, 389.0, 1000.0).  The density and pressure are set          *
 * accordingly.  The gas is assumed to be air.                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState ForwardFacingStep(const HighTemp2D_pState &Wo,
					 const Vector2D &X,
					 const double &h,
					 const double &ho,
					 const double &Re,
					 const double &M) {
  assert(fabs(Re-100.0) < TOLER || fabs(Re-389.0) < TOLER || fabs(Re-1000.0) < TOLER);
  assert(fabs(M-0.20) < TOLER);
  HighTemp2D_pState W;
  int n;
  double eta, c1, c2, c3, c4, c5, f, df, T, To, T1, T2, f1, f2;
  // Set the density.  Quite low eh.
  if (fabs(Re-100.0) < TOLER) W.rho = 0.0046;
  else if (fabs(Re-389.0) < TOLER) W.rho = 0.0104;
  else if (fabs(Re-1000.0) < TOLER) W.rho = 0.045;
  // Determine the temperature.
  c1 = AIR_c1; c2 = AIR_c2; c3 = AIR_c3; c4 = AIR_c4;
  c5 = TWO*W.rho*h*sqrt(W.g*W.R)*M/Re;
  T = 335.0; To = 336.0; n = 0;
  while (fabs(T-To) > 1.0e-06) {
    To = T;
    f = c1*To - c5*(c2 + pow(To,c3) + c4/To);
    df = c1 - c5*(c3*pow(To,c3-ONE) - c4/(To*To));
    T = To - f/df;
    n++;
    if (n > 25) {
      cout << n << " " << T << " " << To << " " << f << " " << df << endl;
      break;
    }
  }
  // Determine the pressure.
  W.p = W.rho*W.R*T;
  // Determine the x-component of the velocity.
  if (X.x > ZERO) eta = ZERO;
  else eta = ONE + TWO*(X.y - ho - h)/h;
  W.v.x = W.a()*M*(ONE - eta*eta);
  // Set the y-component of the velocity.
  W.v.y = ZERO;
  // Return the solution state.
  return W;
}

/**********************************************************************
 * Routine: BC_Characteristic (Characteristic-Based Boundary          *
 *                             Condition)                             *
 *                                                                    *
 * This function returns the boundary solution state for a given      *
 * direction given the primitive solution state on the interior of    *
 * the boundary, Wi, the desired flow state to be imposed at the      *
 * boundary, Wo, and the unit normal vector in the direction of       *
 * interest. The boundary flow state is based on the solution of a    *
 * simplified Riemann problem as described by Gottlieb and Groth      *
 * (1999).  The imposition of the boundary-data respects the          *
 * directions of propogation for the solution characteristics at the  *
 * boundary and thereby ensures that the boundary conditions are not  *
 * ill-posed (i.e., the boundary data is not under- or over-          *
 * determined).                                                       *
 *                                                                    *
 **********************************************************************/

HighTemp2D_pState BC_Characteristic(const HighTemp2D_pState &Wi,
                                 const HighTemp2D_pState &Wo,
	      	                 const Vector2D &norm_dir) {

  HighTemp2D_pState Wi_rotated, Wo_rotated, We, We_rotated;
  char pattern;
  double mi, poi, pai, pbi, mab, mac1, mac2, mbc1, mc1c2, mbd, mad;
  double de, ue, ve, pe, ae, ue_rotated, ve_rotated;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by
  // the unit normal vector.
  Wi_rotated.Rotate(Wi,norm_dir);
  Wo_rotated.Rotate(Wo,norm_dir);

  // Determine the Mach number at the interior node and the imposed
  // boundary to interior node pressure ratio in the rotated frame.
  mi = Wi_rotated.v.x/Wi_rotated.a();
  poi = Wo_rotated.p/Wi_rotated.p;

  // Identify the wave pattern between the interior and boundary node 
  // in the rotated frame.
  if (poi >= ONE) {
    pai = pow(HALF*(Wi_rotated.g+ONE),TWO*Wi_rotated.g*Wi_rotated.gm1i);
    pbi = pow(HALF*(Wi_rotated.g+ONE)+
	      HALF*Wi_rotated.gm1*(Wo_rotated.a()/Wi_rotated.a()), 
	      TWO*Wi_rotated.g*Wi_rotated.gm1i);
    if (poi <= pai) {
      mab = ONE;
      // Pattern A, pi <= po <= pa.
      if (mi >= mab) {
	pattern='a';
	goto impose_boundary_conditions;
      }
      mbc1 = TWO*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
      // Pattern B, pi <= po <= pa.
      if (mi >= mbc1) {
	pattern='b';
	goto impose_boundary_conditions;
      }
      mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
      // Pattern C1, pi <= po <= pa.
      if (mi >= mc1c2) {
	pattern='c';
	goto impose_boundary_conditions;
	// Pattern C2, pi <= po <= pa.
      } else {
	pattern='c';
	goto impose_boundary_conditions;
      }
    } else if (poi <= pbi) {
      mac1 = ONE;
      // Pattern A, pa < po <= pb.
      if (mi >= mac1) {
	pattern='a';
	goto impose_boundary_conditions;
      }
      mbc1 = TWO*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
      mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
      // Pattern C1, pa < po <= pb.
      if (mi >= mc1c2) {
	pattern='c';
	goto impose_boundary_conditions;
	// Pattern C2, pa < po <= pb.
      } else {
	pattern='c';
	goto impose_boundary_conditions;
      }
    } else {
      mac2 = ONE;
      // Pattern A, po > pb.
      if (mi >= mac2) {
	pattern='a';
	goto impose_boundary_conditions;
	// Pattern C2, po > pb.
      } else {
	pattern='c';
	goto impose_boundary_conditions;
      }
    }
  } else {
    mad = ONE;
    // Pattern A, po < pi.
    if (mi >= mad) {
      pattern='a';
      goto impose_boundary_conditions;
    }
    mbd = (Wi_rotated.g+ONE)*Wi_rotated.gm1i*pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-
      TWO*Wi_rotated.gm1i;
    // Pattern D, po < pi.
    if (mi >= mbd) {
      pattern='d';
      goto impose_boundary_conditions;
    }
    mbc1 = TWO*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    // Pattern B, po < pi.
    if (mi >= mbc1) {
      pattern='b';
      goto impose_boundary_conditions;
    }
    mc1c2 = mbc1 - Wo_rotated.a()/Wi_rotated.a();
    // Pattern C1, po < pi.
    if (mi >= mc1c2) {
      pattern='c';
      goto impose_boundary_conditions;
      // Pattern C2, po < pi.
    } else {
      pattern='c';
      goto impose_boundary_conditions;
    }
  }

  // Prescribe the appropriate boundary conditions in the rotated frame.  
 impose_boundary_conditions: ;
  switch(pattern) {
  case 'a' :
    We_rotated.rho   = Wi_rotated.rho;
    We_rotated.v.x   = Wi_rotated.v.x;
    We_rotated.v.y   = Wi_rotated.v.y;
    We_rotated.p     = Wi_rotated.p;
    break;        
  case 'b' :
    We_rotated.rho   = Wi_rotated.rho*pow(poi,ONE/Wi_rotated.g);
    We_rotated.v.x   = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y   = Wi_rotated.v.y;
    We_rotated.p     = Wo_rotated.p;
    break;
  case 'c' :
    We_rotated.rho   = Wo_rotated.rho;
    We_rotated.v.x   = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y   = Wo_rotated.v.y;
    We_rotated.p     = Wo_rotated.p;
    break;
  case 'd' :
    We_rotated.p     = Wi_rotated.p*pow((TWO/(Wi_rotated.g+ONE)+(Wi_rotated.gm1/(Wi_rotated.g+ONE))*mi),TWO*Wi_rotated.g*Wi_rotated.gm1i);
    We_rotated.v.x   = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(We_rotated.p/Wi_rotated.p,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y   = Wi_rotated.v.y;
    We_rotated.rho   = Wi_rotated.rho*pow(We_rotated.p/Wi_rotated.p,ONE/Wi_rotated.g);
    break;
  default:
    We_rotated.rho   = Wo_rotated.rho;
    We_rotated.v.x   = Wo_rotated.v.x;
    We_rotated.v.y   = Wo_rotated.v.y;
    We_rotated.p     = Wo_rotated.p;
    break;
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  We.Rotate(We_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return boundary solution state.
  return We;

}


/**********************************************************************
 * Routine: BC_Characteristic_Pressure                                *
 *   (Characteristic-Based Boundary Condition with Static Pressure    *
 *    Specified Whenever Possible)                                    *
 *                                                                    *
 * This function returns the boundary solution state for a given      *
 * direction given the primitive solution state on the interior of    *
 * the boundary, Wi, the desired flow state to be imposed at the      *
 * boundary, Wo, and the unit normal vector in the direction of       *
 * interest. A simplified characteristic analysis is used to specify  *
 * the boundary flow state in which the static pressure is specified  *
 * whenever possible.  The imposition of the boundary-data respects   *
 * the directions of propogation for the solution characteristics at  *
 * the boundary and thereby ensures that the boundary conditions are  *
 * not ill-posed (i.e., the boundary data is not under- or over-      *
 * determined).  The following procedure is adopted:                  *
 *                                                                    *
 * 1) for supersonic outflow: constant extrapolation is employed to   *
 *    specify boundary conditions using known interior solution       *
 *    state,                                                          *
 * 2) for subsonic outflow: boundary conditions specified by          *
 *    employing a 1D unsteady isentropic wave approximation to match  *
 *    boundary pressure,                                              *
 * 3) for subsonic inflow: boundary conditions specified by employing *
 *    a 1D unsteady isentropic wave approximation to match both       *
 *    boundary state pressure and sound speed,                        *
 * 4) for supersonic inflow: the known boundary state is used to      *
 *    specify the boundary state.                                     *
 *                                                                    *
 **********************************************************************/

HighTemp2D_pState BC_Characteristic_Pressure(const HighTemp2D_pState &Wi,
						 const HighTemp2D_pState &Wo,
						 const Vector2D &norm_dir) {

  HighTemp2D_pState Wi_rotated, Wo_rotated, Wb, Wb_rotated;
  double mi, ab;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by 
  // the unit normal vector.
  Wi_rotated.Rotate(Wi,norm_dir);
  Wo_rotated.Rotate(Wo,norm_dir);

  // Determine the Mach number at the interior node.
  mi = Wi_rotated.v.x/Wi_rotated.a();

  // Boundary condition for supersonic outflow.
  if (mi >= ONE) {
    Wb_rotated.rho   = Wi_rotated.rho;
    Wb_rotated.v.x   = Wi_rotated.v.x;
    Wb_rotated.v.y   = Wi_rotated.v.y;
    Wb_rotated.p     = Wi_rotated.p;
    Wb_rotated.k     = Wi_rotated.k;
    Wb_rotated.omega = Wi_rotated.omega;

  // Boundary condition for subsonic outflow.  Pressure specified.
  } else if (mi >= ZERO) {
    Wb_rotated.p   = Wo_rotated.p;
    Wb_rotated.rho = Wi_rotated.rho*pow(Wb_rotated.p/Wi_rotated.p,ONE/Wi_rotated.g);
    ab = sqrt(Wi_rotated.g*Wb_rotated.p/Wb_rotated.rho);
    Wb_rotated.v.x = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wi_rotated.gm1i;
    Wb_rotated.v.y = Wi_rotated.v.y;
    Wb_rotated.k     = Wo_rotated.k;
    Wb_rotated.omega = Wo_rotated.omega;

  // Boundary condition for subsonic inflow.  Pressure specified.
  } else if (mi >= -ONE) {
    Wb_rotated.p   = Wo_rotated.p;
    Wb_rotated.rho = Wo_rotated.rho;
    ab = sqrt(Wo_rotated.g*Wb_rotated.p/Wb_rotated.rho);
    Wb_rotated.v.x = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wo_rotated.gm1i;
    Wb_rotated.v.y = Wo_rotated.v.y;
    Wb_rotated.k     = Wo_rotated.k;
    Wb_rotated.omega = Wo_rotated.omega;

  // Boundary condition for supersonic inflow.
  } else {
    Wb_rotated.rho = Wo_rotated.rho;
    Wb_rotated.v.x = Wo_rotated.v.x;
    Wb_rotated.v.y = Wo_rotated.v.y;
    Wb_rotated.p   = Wo_rotated.p;
    Wb_rotated.k     = Wo_rotated.k;
    Wb_rotated.omega = Wo_rotated.omega;
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  Wb.Rotate(Wb_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return boundary solution state.
  return Wb;

}

/**********************************************************************
 * Routine: BC_Characteristic_Mach_Number                             *
 *   (Characteristic-Based Boundary Condition with Flow Mach Number   *
 *    Specified Whenever Possible)                                    *
 *                                                                    *
 * This function returns the boundary solution state for a given      *
 * direction given the primitive solution state on the interior of    *
 * the boundary, Wi, the desired flow state to be imposed at the      *
 * boundary, Wo, and the unit normal vector in the direction of       *
 * interest. A simplified characteristic analysis is used to specify  *
 * the boundary flow state in which the flow Mach number is specified *
 * whenever possible.                                                 *
 * The imposition of the boundary-data respects the directions of     *
 * propogation for the solution characteristics at the boundary and   *
 * thereby ensures that the boundary conditions are not ill-posed     *
 * (i.e., the boundary data is not under- or over-determined).  The   *
 * following procedure is adopted:                                    *
 *                                                                    *
 * 1) for supersonic outflow: constant extrapolation is employed to   *
 *    specify boundary conditions using known interior solution       *
 *    state,                                                          *
 * 2) for subsonic outflow: boundary conditions specified by          *
 *    employing a 1D unsteady isentropic wave approximation to match  *
 *    boundary Mach number,                                           *
 * 3) for subsonic inflow: boundary conditions specified by employing *
 *    a 1D unsteady isentropic wave approximation to match both       *
 *    boundary state Mach number and density,                         *
 * 4) for supersonic inflow: the known boundary state is used to      *
 *    specify the boundary state.                                     *
 *                                                                    *
 **********************************************************************/


HighTemp2D_pState BC_Characteristic_Mach_Number(const HighTemp2D_pState &Wi,
                                             const HighTemp2D_pState &Wo,
	      	                             const Vector2D &norm_dir) {

  HighTemp2D_pState Wi_rotated, Wo_rotated, Wb, Wb_rotated;
  double mi, mb, ab;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by 
  // the unit normal vector.
  Wi_rotated.Rotate(Wi,norm_dir);
  Wo_rotated.Rotate(Wo,norm_dir);

  // Determine the Mach number at the interior node.
  mi = Wi_rotated.v.x/Wi_rotated.a();

  // Boundary condition for supersonic outflow.
  if (mi >= ONE) {
    Wb_rotated.rho = Wi_rotated.rho;
    Wb_rotated.v.x = Wi_rotated.v.x;
    Wb_rotated.v.y = Wi_rotated.v.y;
    Wb_rotated.p   = Wi_rotated.p;

  // Boundary condition for subsonic outflow. Mach number specified.
  } else if (mi >= ZERO) {
    Wb_rotated.rho = Wi_rotated.rho;
    mb = Wo_rotated.v.x/Wo_rotated.a();
    ab = (Wi_rotated.v.x+TWO*Wi_rotated.a()*Wi_rotated.gm1i)/(mb+TWO*Wi_rotated.gm1i);
    Wb_rotated.p   = Wb_rotated.rho*ab*ab/Wi_rotated.g;
    Wb_rotated.v.x = mb*ab;
    Wb_rotated.v.y = Wi_rotated.v.y;

  // Boundary condition for subsonic inflow.  Pressure specified.
  } else if (mi >= -ONE) {
    Wb_rotated.rho = Wo_rotated.rho;
    mb = Wo_rotated.v.x/Wo_rotated.a();
    ab = (Wi_rotated.v.x+TWO*Wi_rotated.a()*Wi_rotated.gm1i)/(mb+TWO*Wi_rotated.gm1i);
    Wb_rotated.p   = Wb_rotated.rho*ab*ab/Wo_rotated.g;
    Wb_rotated.v.x = mb*ab;
    Wb_rotated.v.y = Wo_rotated.v.y;    

  // Boundary condition for supersonic inflow.
  } else {
    Wb_rotated.rho = Wo_rotated.rho;
    Wb_rotated.v.x = Wo_rotated.v.x;
    Wb_rotated.v.y = Wo_rotated.v.y;
    Wb_rotated.p   = Wo_rotated.p;
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  Wb.Rotate(Wb_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return boundary solution state.
  return Wb;

}


/**********************************************************************
 * Routine: WaveSpeedPos                                              *
 *                                                                    *
 * This function returns the positive parts of the elemental wave     *
 * speeds or eigenvalues.                                             *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState WaveSpeedPos(const HighTemp2D_pState &lambdas_a,
				   const HighTemp2D_pState &lambdas_l,
				   const HighTemp2D_pState &lambdas_r) {
  HighTemp2D_pState Wsp;
  for (int n = 0; n < NUM_VAR_HIGHTEMP2D; n++)
    Wsp[n] = HALF*(lambdas_a[n]+fabs(lambdas_a[n]));
  return Wsp;
}

/**********************************************************************
 * Routine: WaveSpeedNeg                                              *
 *                                                                    *
 * This function returns the negative parts of the elemental wave     *
 * speeds or eigenvalues.                                             *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState WaveSpeedNeg(const HighTemp2D_pState &lambdas_a,
                            const HighTemp2D_pState &lambdas_l,
                            const HighTemp2D_pState &lambdas_r) {
  HighTemp2D_pState Wsn;
  for (int n = 0; n < NUM_VAR_HIGHTEMP2D; n++)
    Wsn[n] = HALF*(lambdas_a[n]-fabs(lambdas_a[n]));
  return Wsn;
}

/**********************************************************************
 * Routine: WaveSpeedAbs                                              *
 *                                                                    *
 * This function returns the absolute values of the elemental wave    *
 * speeds or eigenvalues.                                             *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState WaveSpeedAbs(const HighTemp2D_pState &lambdas_a,
				   const HighTemp2D_pState &lambdas_l,
				   const HighTemp2D_pState &lambdas_r) {
  return HighTemp2D_pState(fabs(lambdas_a[1]),
			       fabs(lambdas_a[2]),
			       fabs(lambdas_a[3]),
			       fabs(lambdas_a[4]),
			       fabs(lambdas_a[5]),
			       fabs(lambdas_a[6]));
}

/**********************************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the positive parts of the corrected          *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState HartenFixPos(const HighTemp2D_pState &lambdas_a,
				   const HighTemp2D_pState &lambdas_l,
				   const HighTemp2D_pState &lambdas_r) {
  return HighTemp2D_pState(HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]),
			       HALF*(lambdas_a[2]+fabs(lambdas_a[2])),
			       HALF*(lambdas_a[3]+fabs(lambdas_a[3])),
			       HartenFixPos(lambdas_a[4],lambdas_l[4],lambdas_r[4]),
			       HALF*(lambdas_a[5]+fabs(lambdas_a[5])),
			       HALF*(lambdas_a[6]+fabs(lambdas_a[6])));
}

/**********************************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the negative parts of the corrected          *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState HartenFixNeg(const HighTemp2D_pState &lambdas_a,
				   const HighTemp2D_pState &lambdas_l,
				   const HighTemp2D_pState &lambdas_r) {
  return HighTemp2D_pState(HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]),
			       HALF*(lambdas_a[2]-fabs(lambdas_a[2])),
			       HALF*(lambdas_a[3]-fabs(lambdas_a[3])),
			       HartenFixNeg(lambdas_a[4],lambdas_l[4],lambdas_r[4]),
			       HALF*(lambdas_a[5]-fabs(lambdas_a[5])),
			       HALF*(lambdas_a[6]-fabs(lambdas_a[6])));
}

/**********************************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the absolute values of the corrected         *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState HartenFixAbs(const HighTemp2D_pState &lambdas_a,
				   const HighTemp2D_pState &lambdas_l,
				   const HighTemp2D_pState &lambdas_r) {
  return HighTemp2D_pState(HartenFixAbs(lambdas_a[1],lambdas_l[1],lambdas_r[1]),
			       fabs(lambdas_a[2]),
			       fabs(lambdas_a[3]),
			       HartenFixAbs(lambdas_a[4],lambdas_l[4],lambdas_r[4]),
			       fabs(lambdas_a[5]),
			       fabs(lambdas_a[6]));
}

/**********************************************************************
 * Routine: FluxGodunov_n (Godunov's flux function, n-direction)      *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate     *
 * the problem to a local frame aligned with the unit normal vector   *
 * and then solving exactly the Riemann problem associated with the   *
 * two states in the rotated frame.  See Godunov (1959).              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGodunov_n(const HighTemp2D_pState &Wl,
				    const HighTemp2D_pState &Wr,
				    const Vector2D &norm_dir) {

  HighTemp2D_pState W, Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;
  
  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);
  
  // Solve the Riemann problem in the rotated frame and evaluate the 
  // intermediate solution state and flux.
  W = Riemann(Wl_rotated,Wr_rotated);
  //THIS IS W.F() in the NS code!
  Flux_rotated = W.Fx();
  
  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

HighTemp2D_cState FluxGodunov_n(const HighTemp2D_cState &Ul,
				    const HighTemp2D_cState &Ur,
				    const Vector2D &norm_dir) {
  return FluxGodunov_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxGodunov_MB_n (Godunov's flux function, n-direction)   *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate     *
 * the problem to a local frame aligned with the unit normal vector   *
 * and then solving exactly the Riemann problem associated with the   *
 * two states in the rotated frame.  See Godunov (1959).              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGodunov_MB_n(const HighTemp2D_pState &Wl,
				       const HighTemp2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir) {

  HighTemp2D_pState W_rotated, Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Solve the Riemann problem in the rotated frame and evaluate the 
  // intermediate state solution flux.
  W_rotated = Riemann(Wl_rotated,Wr_rotated);

  // Transform back into the moving body frame of reference.
  W_rotated.v.x += V_rotated.x;

  // Evaluate the intermediate state solution flux in the rotated frame.
  //THIS IS .F(V_rotated) in NS2D code
  Flux_rotated = W_rotated.Fx(V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: StateGodunov_n (Godunov's flux function, n-direction)     *
 *                                                                    *
 * This function returns the intermediate solution state for an       *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate     *
 * the problem to a local frame aligned with the unit normal vector   *
 * and then solving exactly the Riemann problem associated with the   *
 * two states in the rotated frame.  See Godunov (1959).              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_pState StateGodunov_n(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr,
				     const Vector2D &norm_dir) {

  HighTemp2D_pState W_rotated, Wl_rotated, Wr_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Solve the Riemann problem in the rotated frame and evaluate the 
  // intermediate solution state.
  W_rotated = Riemann(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution state.
	HighTemp2D_pState Wtmp; Wtmp.Rotate(W_rotated,Vector2D(norm_dir.x,-norm_dir.y));
	return Wtmp;

}

/**********************************************************************
 * Routine: FluxRoe (Roe's flux function)                             *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * given left and right solution states by using the "linearized"     *
 * approximate Riemann solver of Roe for the two states.  See Roe     *
 * (1981).                                                            *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRoe(const HighTemp2D_pState &Wl,
			      const HighTemp2D_pState &Wr) {

  HighTemp2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_HIGHTEMP2D : NUM_VAR_HIGHTEMP2D-2;
 
  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the primitive solution states.
  dWrl = Wr - Wl;

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();
  
  // Determine the intermediate state flux.
  if (Wa.v.x >= ZERO) {
    Flux = Wl.Fx();
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] < ZERO) {
	Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
      }
    }
  } else {
    Flux = Wr.Fx();
    wavespeeds = HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] > ZERO) {
	Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
      }
    }
  }

  // Return the solution flux.
   return Flux;

}

HighTemp2D_cState FluxRoe(const HighTemp2D_cState &Ul,
			      const HighTemp2D_cState &Ur) {
  return FluxRoe(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxRoe_n (Roe's flux function, n-direction)              *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the "linearized" approximate Riemann solver of Roe   *
 * to specify the flux in terms of the rotated solution states.  See  *
 * Roe (1981).                                                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRoe_n(const HighTemp2D_pState &Wl,
				const HighTemp2D_pState &Wr,
				const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRoe(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxRoe_n(const HighTemp2D_cState &Ul,
				const HighTemp2D_cState &Ur,
				const Vector2D &norm_dir) {
  return FluxRoe_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxGlaister (modiefied Roe's flux function)              *
 *  FOR HIGH-TEMPERATURE AIR EFFECTS                                  *
 * This function returns the intermediate state solution flux for the *
 * given left and right solution states by using the "linearized"     *
 * approximate Riemann solver of Roe for the two states.  See Roe     *
 * (1981) and Glaister (1987).                                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGlaister(const HighTemp2D_pState &Wl,
			      const HighTemp2D_pState &Wr) {

  HighTemp2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_HIGHTEMP2D : NUM_VAR_HIGHTEMP2D-2;

  //Glaister-algorithm definitions
  int deltaE, deltaRHO;
  double pspl, psml, pspr, psmr, prl, plr, plsp, plsm, prsp, prsm;
  double de, drho, dpde, dpdrho, dpder, dpdel, dpdrhol, dpdrhor;
  double RHOstarP, RHOstarM, EstarP, EstarM;

  // Evaluate the modified Roe-average (Glaister) primitive solution state.
  //Wa = RoeAverage(Wl,Wr);
  double hl, hr, el, er, srhol, srhor, ha, ea, cAvg, cl_local, cr_local;
  Wa.rho = ZERO; 
  Wa.p = ZERO; //Wa.Vacuum();

  // Determine the left and right state specific enthalpies and square
  // roots of the density.
  el    = Wl.e();
  er    = Wr.e();
  hl    = el + Wl.p/Wl.rho;
  hr    = er + Wr.p/Wr.rho;
  //el    = Wl.e();
  // er    = Wr.e();
  srhol = sqrt(Wl.rho);
  srhor = sqrt(Wr.rho);
  
  // Determine the appropriate Roe averages.
  Wa.rho   = srhol*srhor;
  Wa.v.x   = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
  Wa.v.y   = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
  Wa.k     = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
  Wa.omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
  ha = (srhol*hl + srhor*hr)/(srhol+srhor);
  // Glaister modification
  ea = (srhol*el + srhor*er)/(srhol+srhor);
  Wa.p = Wa.rho*(ha - ea);
  //Wa.p = Tgas_p(ea,Wa.rho);  

  //calculation of modified sound speed for Average State  
  RHOstarP = Wa.rho*HTONEPLUST;
  RHOstarM = Wa.rho*HTONEMINT;
  EstarP = ea*HTONEPLUST;
  EstarM = ea*HTONEMINT;
  de = (er-el);
  drho = (Wr.rho-Wl.rho);
 
  if (abs(de)<(TOLER3*ea))
    deltaE = 0;
  else deltaE = 1;
  
  if (abs(drho)<(TOLER3*Wa.rho))
    deltaRHO = 0;
  else deltaRHO = 1;

  //Find dpde derivative
  if (deltaE == 1){
    prl = Tgas_p(er, Wl.rho);
    plr = Tgas_p(el, Wr.rho);
    dpde = 0.5*(Wr.p - plr + prl - Wl.p)/de;
  }else{
    pspl  = Tgas_p(EstarP,Wl.rho);
    psml  = Tgas_p(EstarM,Wl.rho);
    dpdel = (pspl-psml)/(2.00*HTTOL*ea);       
    pspr  = Tgas_p(EstarP,Wr.rho);
    psmr  = Tgas_p(EstarM,Wr.rho);
    dpder = (pspr-psmr)/(2.00*HTTOL*ea);    
    dpde  = 0.50*(dpdel+dpder);
  }

  //Find dpdrho derivative
   if (deltaRHO == 1){
     if (deltaE == 1)
       dpdrho = 0.5*(Wr.p - prl + plr - Wl.p)/drho;
     else {
       prl = Tgas_p(er, Wl.rho);
       plr = Tgas_p(el, Wr.rho);
       dpdrho = 0.5*(Wr.p - prl + plr - Wl.p)/drho;
     }
   }else{
    plsp    = Tgas_p(el,RHOstarP);
    plsm    = Tgas_p(el,RHOstarM); 
    dpdrhol = (plsp-plsm)/(2.00*HTTOL*Wa.rho);
    prsp    = Tgas_p(er,RHOstarP);
    prsm    = Tgas_p(er,RHOstarM);    
    dpdrhor = (prsp-prsm)/(2.00*HTTOL*Wa.rho);
    dpdrho  = 0.50*(dpdrhol + dpdrhor);
  }
  //Glaister approximation for Average MODIFIED sound-speed 
  //aAvg = sqrt(Wa.p*dpde/(Wa.rho*Wa.rho)+dpdrho); 
   double aa2;
  aa2 = Wa.p*dpde/(Wa.rho*Wa.rho)+dpdrho + 2.0*Wa.k/3.0 + (2.0*Wa.k*dpde)/(3.0*Wa.rho);   
  cAvg = sqrt(Wa.p*dpde/(Wa.rho*Wa.rho)+dpdrho + 2.0*Wa.k/3.0 + (2.0*Wa.k*dpde)/(3.0*Wa.rho));   

 //added from IDEAL - just for testing
  // Wa.gm1 = GAMMA_AIR-ONE;  
  //Wa.g = GAMMA_AIR;
  //aa2 = Wa.gm1*(ha - HALF*(sqr(Wa.v.x) + sqr(Wa.v.y)) - Wa.k);
  //Wa.p = Wa.rho*aa2/Wa.g;
  //cAvg = sqrt(aa2);

  // Evaluate the jumps in the primitive solution states.
  dWrl = Wr - Wl;

  // Evaluate the left, right, and average state eigenvalues.
  cl_local = Wl.c();
  cr_local = Wr.c();
  lambdas_l = Wl.lambda_x(cl_local);
  lambdas_r = Wr.lambda_x(cr_local);
  //using Glaister's computed sound-speed
  lambdas_a = Wa.lambda_x(cAvg);
 if (Wa.v.x >= ZERO) {
    Flux = Wl.Fx();
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] < ZERO) {
	Flux += wavespeeds[i]*(Wa.lp_x(i,cAvg)*dWrl)*Wa.rc_x(i,dpde,dpdrho,cAvg);
      }
    }
  } else {
    Flux = Wr.Fx();
    wavespeeds = HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] > ZERO) {
       	Flux -= wavespeeds[i]*(Wa.lp_x(i,cAvg)*dWrl)*Wa.rc_x(i,dpde,dpdrho,cAvg);
      }
    }
  }

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxGlaister(const HighTemp2D_cState &Ul,
			      const HighTemp2D_cState &Ur) {
  return FluxGlaister(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: Flux_n (Roe's flux function, n-direction)              *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the "linearized" approximate Riemann solver of Roe   *
 * to specify the flux in terms of the rotated solution states.  See  *
 * Roe (1981).                                                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGlaister_n(const HighTemp2D_pState &Wl,
				const HighTemp2D_pState &Wr,
				const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxGlaister(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxGlaister_n(const HighTemp2D_cState &Ul,
				const HighTemp2D_cState &Ur,
				const Vector2D &norm_dir) {
  return FluxGlaister_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxRoe_MB (Roe's flux function)                          *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * given left and right solution states by using the "linearized"     *
 * approximate Riemann solver of Roe for the two states.  See Roe     *
 * (1981).                                                            *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRoe_MB(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &V) {

  HighTemp2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_HIGHTEMP2D : NUM_VAR_HIGHTEMP2D-2;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the primitive solution states.
  dWrl = Wr - Wl;

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x(V);
  lambdas_r = Wr.lambda_x(V);
  lambdas_a = Wa.lambda_x(V);

  // Determine the intermediate state flux.
  if (Wa.v.x >= ZERO) {
    Flux = Wl.Fx(V);
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] < ZERO)
	Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  } else {
    Flux = Wr.Fx(V);
    wavespeeds = HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] > ZERO)
	Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  }

  // Return the solution flux.
  return Flux;
    
}

HighTemp2D_cState FluxRoe_MB(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &V) {
  return FluxRoe_MB(Ul.W(),Ur.W(),V);
}

/**********************************************************************
 * Routine: FluxRoe_MB_n (Roe's flux function, n-direction)           *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the "linearized" approximate Riemann solver of Roe   *
 * to specify the flux in terms of the rotated solution states.  See  *
 * Roe (1981).                                                        *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRoe_MB_n(const HighTemp2D_pState &Wl,
				   const HighTemp2D_pState &Wr,
				   const Vector2D &V,
				   const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_pState Wl_nonstationary, Wr_nonstationary;
  HighTemp2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Move left and right states back into nonstationary frame.
  Wl_nonstationary = Translate(Wl,-V);
  Wr_nonstationary = Translate(Wr,-V);

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl_nonstationary,norm_dir);
  Wr_rotated.Rotate(Wr_nonstationary,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRoe_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxRoe_MB_n(const HighTemp2D_cState &Ul,
				   const HighTemp2D_cState &Ur,
				   const Vector2D &V,
				   const Vector2D &norm_dir) {
  return FluxRoe_MB_n(Ul.W(),Ur.W(),V,norm_dir);
}

/**********************************************************************
 * Routine: FluxRusanov (Rusanov flux function, x-direction)          *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * Rusanov approximation for the fluxes.  See Rusanov (1964).         *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRusanov(const HighTemp2D_pState &Wl,
				  const HighTemp2D_pState &Wr) {

  double wavespeed_max;
  HighTemp2D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux, dUrl;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_HIGHTEMP2D : NUM_VAR_HIGHTEMP2D-2;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U()-Wl.U();

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();

  // Determine the intermediate state flux.
  Flux = HALF*(Wl.Fx()+Wr.Fx());
  wavespeeds = HartenFixAbs(lambdas_a,lambdas_l,lambdas_r);

  wavespeed_max = wavespeeds[1];
  for (int i = 2; i <= NumVar; i++)
    wavespeed_max = max(wavespeed_max,wavespeeds[i]);
  
  Flux -= HALF*wavespeed_max*dUrl;

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxRusanov(const HighTemp2D_cState &Ul,
				  const HighTemp2D_cState &Ur) {
  return FluxRusanov(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxRusanov_n (Rusanov flux function, n-direction)        *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the Rusanov approximation for the intermediate state *
 * flux in terms of the rotated solution states.  See Rusanov (1964). *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxRusanov_n(const HighTemp2D_pState &Wl,
				    const HighTemp2D_pState &Wr,
				    const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRusanov(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxRusanov_n(const HighTemp2D_cState &Ul,
				    const HighTemp2D_cState &Ur,
				    const Vector2D &norm_dir) {
  return FluxRusanov_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxHLLE (Harten-Lax-van Leer flux function, x-direction) *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * so-called Harten-Lax-van Leer approximation for the fluxes.  See   *
 * Harten, Lax, van Leer (1983).                                      *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLE(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  double wavespeed_l, wavespeed_r;
  HighTemp2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux, dUrl;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U() - Wl.U();

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();
  
  // Determine the intermediate state flux.
  wavespeed_l = min(lambdas_l[1],lambdas_a[1]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
  
  wavespeed_l = min(wavespeed_l,ZERO);
  wavespeed_r = max(wavespeed_r,ZERO);

  if (wavespeed_l >= ZERO) {
    Flux = Wl.Fx();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.Fx();
  } else {
    Flux = (((wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
	     + (wavespeed_l*wavespeed_r)*dUrl)/
	    (wavespeed_r-wavespeed_l));
  }

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxHLLE(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxHLLE(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux function,            *
 *                      n-direction)                                  *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the so-called Harten-Lax-van Leer approximation to   *
 * specify the intermediate state fluxes in terms of the rotated      *
 * solution states.  See Harten, Lax, van Leer (1983).                *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLE_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;
  
  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated
  // frame.
  Flux_rotated = FluxHLLE(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and
  // return the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxHLLE_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxHLLE_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxGHLLE (Glaister, Harten-Lax-van Leer flux function,   *
 *                                                       x-direction) *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * so-called Harten-Lax-van Leer approximation for the fluxes.  See   *
 * Harten, Lax, van Leer (1983). The average state is computed using  *
 * Glaister's modification for high temperature (1988).               *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGHLLE(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  double wavespeed_l, wavespeed_r;
  HighTemp2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux, dUrl;

 //Glaister-algorithm definitions
  int deltaE, deltaRHO;
  double pspl, psml, pspr, psmr, prl, plr, plsp, plsm, prsp, prsm;
  double de, drho, dpde, dpdrho, dpder, dpdel, dpdrhol, dpdrhor;
  double RHOstarP, RHOstarM, EstarP, EstarM;
  double hl, hr, el, er, srhol, srhor, ha, ea, cAvg, cl_local, cr_local;
  Wa.rho = ZERO; 
  Wa.p = ZERO; //Wa.Vacuum();
 
  // Evaluate the modified Roe-average (Glaister) primitive solution state.

  // Determine the left and right state specific enthalpies and square
  // roots of the density.
  el    = Wl.e();
  er    = Wr.e();
  hl    = el + Wl.p/Wl.rho;
  hr    = er + Wr.p/Wr.rho;
  srhol = sqrt(Wl.rho);
  srhor = sqrt(Wr.rho);
  
  // Determine the appropriate Roe averages.
  Wa.rho   = srhol*srhor;
  Wa.v.x   = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
  Wa.v.y   = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
  Wa.k     = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
  Wa.omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
  ha = (srhol*hl + srhor*hr)/(srhol+srhor);
  // Glaister modification
  ea = (srhol*el + srhor*er)/(srhol+srhor);
  Wa.p = Wa.rho*(ha - ea);
 
  //calculation of modified sound speed for Average State  
  RHOstarP = Wa.rho*HTONEPLUST;
  RHOstarM = Wa.rho*HTONEMINT;
  EstarP = ea*HTONEPLUST;
  EstarM = ea*HTONEMINT;
  de = (er-el);
  drho = (Wr.rho-Wl.rho);
 
  if (abs(de)<(TOLER3*ea))
    deltaE = 0;
  else deltaE = 1;
  
  if (abs(drho)<(TOLER3*Wa.rho))
    deltaRHO = 0;
  else deltaRHO = 1;

  //Find dpde derivative
  if (deltaE == 1){
    prl = Tgas_p(er, Wl.rho);
    plr = Tgas_p(el, Wr.rho);
    dpde = 0.5*(Wr.p - plr + prl - Wl.p)/de;
  }else{
    pspl  = Tgas_p(EstarP,Wl.rho);
    psml  = Tgas_p(EstarM,Wl.rho);
    dpdel = (pspl-psml)/(2.00*HTTOL*ea);       
    pspr  = Tgas_p(EstarP,Wr.rho);
    psmr  = Tgas_p(EstarM,Wr.rho);
    dpder = (pspr-psmr)/(2.00*HTTOL*ea);    
    dpde  = 0.50*(dpdel+dpder);
  }

  //Find dpdrho derivative
   if (deltaRHO == 1){
      if (deltaE == 1)
       dpdrho = 0.5*(Wr.p - prl + plr - Wl.p)/drho;
     else {
       prl = Tgas_p(er, Wl.rho);
       plr = Tgas_p(el, Wr.rho);
       dpdrho = 0.5*(Wr.p - prl + plr - Wl.p)/drho;
     }
  }else{
    plsp    = Tgas_p(el,RHOstarP);
    plsm    = Tgas_p(el,RHOstarM); 
    dpdrhol = (plsp-plsm)/(2.00*HTTOL*Wa.rho);
    prsp    = Tgas_p(er,RHOstarP);
    prsm    = Tgas_p(er,RHOstarM);    
    dpdrhor = (prsp-prsm)/(2.00*HTTOL*Wa.rho);
    dpdrho  = 0.50*(dpdrhol + dpdrhor);
  }
  cAvg = sqrt(Wa.p*dpde/(Wa.rho*Wa.rho)+dpdrho + 2.0*Wa.k/3.0 + (2.0*Wa.k*dpde)/(3.0*Wa.rho));   

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U() - Wl.U();

  // Evaluate the left, right, and average state eigenvalues.
  cl_local = Wl.c();
  cr_local = Wr.c();
  lambdas_l = Wl.lambda_x(cl_local);
  lambdas_r = Wr.lambda_x(cr_local);
  //using Glaister's computed sound-speed
  lambdas_a = Wa.lambda_x(cAvg);
  
  // Determine the intermediate state flux.
  wavespeed_l = min(lambdas_l[1],lambdas_a[1]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
  
  wavespeed_l = min(wavespeed_l,ZERO);
  wavespeed_r = max(wavespeed_r,ZERO);

  //cout<<"Computing inside GHLLE function :) "<<endl; 

  if (wavespeed_l >= ZERO) {
    Flux = Wl.Fx();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.Fx();
  } else {
    Flux = (((wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
	     + (wavespeed_l*wavespeed_r)*dUrl)/
	    (wavespeed_r-wavespeed_l));
  }

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxGHLLE(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxGHLLE(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux function,            *
 *                      n-direction)                                  *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the so-called Harten-Lax-van Leer approximation to   *
 * specify the intermediate state fluxes in terms of the rotated      *
 * solution states.  See Harten, Lax, van Leer (1983).                *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGHLLE_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;
  
  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated
  // frame.
  Flux_rotated = FluxGHLLE(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and
  // return the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxGHLLE_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxGHLLE_n(Ul.W(),Ur.W(),norm_dir);
}

/*********************************************************
 * Routine: HLLE_wavespeeds                              *
 *                                                       *
 * This function returns lambda plus and lambda minus    *
 * for rotated Riemann problem aligned with norm_dir     *
 * given unrotated solution states Wl and Wr.             *
 * Note: wavespeed.x = wavespeed_l = lambda minus.       *
 *       wavespeed.y = wavespeed_r = lambda plus.        *
 *                                                       *
 *********************************************************/
Vector2D HLLE_wavespeeds(
		const HighTemp2D_pState &Wl,
		const HighTemp2D_pState &Wr,
		const Vector2D &norm_dir) 
{
    Vector2D wavespeed;
    HighTemp2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  

    /* Use rotated values to calculate eignvalues */
    Wl_n.Rotate(Wl, norm_dir);
    Wr_n.Rotate(Wr, norm_dir);

    /* Evaluate the Roe-average primitive solution state. */
    Wa_n = RoeAverage(Wl_n, Wr_n);
    
    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl_n.lambda_x();
    lambdas_r = Wr_n.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed.x = min(lambdas_l[1], lambdas_a[1]);
    wavespeed.y = max(lambdas_r[4], lambdas_a[4]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    return (wavespeed);

}

/*********************************************************
 * Routine: GHLLE_wavespeeds                              *
 *                                                       *
 * This function returns lambda plus and lambda minus    *
 * for rotated Riemann problem aligned with norm_dir     *
 * given unrotated solution states Wl and Wr.             *
 * Note: wavespeed.x = wavespeed_l = lambda minus.       *
 *       wavespeed.y = wavespeed_r = lambda plus.        *
 *                                                       *
 *********************************************************/
Vector2D GHLLE_wavespeeds(
		const HighTemp2D_pState &Wl,
		const HighTemp2D_pState &Wr,
		const Vector2D &norm_dir) 
{
    Vector2D wavespeed;
    HighTemp2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  
     //Glaister-algorithm definitions
    int deltaE, deltaRHO;
    double pspl, psml, pspr, psmr, prl, plr, plsp, plsm, prsp, prsm;
    double de, drho, dpde, dpdrho, dpder, dpdel, dpdrhol, dpdrhor;
    double RHOstarP, RHOstarM, EstarP, EstarM;
    double hl, hr, el, er, srhol, srhor, ha, ea, cAvg;

    /* Use rotated values to calculate eignvalues */
    Wl_n.Rotate(Wl, norm_dir);
    Wr_n.Rotate(Wr, norm_dir);

    /* Evaluate the Roe-average primitive solution state. */
    //Wa_n = RoeAverage(Wl_n, Wr_n);
  
    Wa_n.rho = ZERO; 
    Wa_n.p = ZERO; //Wa.Vacuum();
 
  // Evaluate the modified Roe-average (Glaister) primitive solution state.

  // Determine the left and right state specific enthalpies and square
  // roots of the density.
  hl    = Wl_n.e() + Wl_n.p/Wl_n.rho;
  hr    = Wr_n.e() + Wr_n.p/Wr_n.rho;
  el    = Wl_n.e();
  er    = Wr_n.e();
  srhol = sqrt(Wl_n.rho);
  srhor = sqrt(Wr_n.rho);
  
  // Determine the appropriate Roe averages.
  Wa_n.rho   = srhol*srhor;
  Wa_n.v.x   = (srhol*Wl_n.v.x + srhor*Wr_n.v.x)/(srhol+srhor);
  Wa_n.v.y   = (srhol*Wl_n.v.y + srhor*Wr_n.v.y)/(srhol+srhor);
  Wa_n.k     = (srhol*Wl_n.k + srhor*Wr_n.k)/(srhol+srhor);
  Wa_n.omega = (srhol*Wl_n.omega + srhor*Wr_n.omega)/(srhol+srhor);
  ha = (srhol*hl + srhor*hr)/(srhol+srhor);
  // Glaister modification
  ea = (srhol*el + srhor*er)/(srhol+srhor);
  Wa_n.p = Wa_n.rho*(ha - ea);
 
  //calculation of modified sound speed for Average State  
  RHOstarP = Wa_n.rho*HTONEPLUST;
  RHOstarM = Wa_n.rho*HTONEMINT;
  EstarP = ea*HTONEPLUST;
  EstarM = ea*HTONEMINT;
  de = (er-el);
  drho = (Wr_n.rho-Wl_n.rho);
 
  if (abs(de)<(TOLER3*ea))
    deltaE = 0;
  else deltaE = 1;
  
  if (abs(drho)<(TOLER3*Wa_n.rho))
    deltaRHO = 0;
  else deltaRHO = 1;

  //Find dpde derivative
  if (deltaE == 1){
    prl = Tgas_p(er, Wl_n.rho);
    plr = Tgas_p(el, Wr_n.rho);
    dpde = 0.5*(Wr_n.p - plr + prl - Wl_n.p)/de;
  }else{
    pspl  = Tgas_p(EstarP,Wl_n.rho);
    psml  = Tgas_p(EstarM,Wl_n.rho);
    dpdel = (pspl-psml)/(2.00*HTTOL*ea);       
    pspr  = Tgas_p(EstarP,Wr_n.rho);
    psmr  = Tgas_p(EstarM,Wr_n.rho);
    dpder = (pspr-psmr)/(2.00*HTTOL*ea);    
    dpde  = 0.50*(dpdel+dpder);
  }

  //Find dpdrho derivative
   if (deltaRHO == 1){
    prl = Tgas_p(er, Wl_n.rho);
    plr = Tgas_p(el, Wr_n.rho);
    dpdrho = 0.5*(Wr_n.p - prl + plr - Wl_n.p)/drho;
  }else{
    plsp    = Tgas_p(el,RHOstarP);
    plsm    = Tgas_p(el,RHOstarM); 
    dpdrhol = (plsp-plsm)/(2.00*HTTOL*Wa_n.rho);
    prsp    = Tgas_p(er,RHOstarP);
    prsm    = Tgas_p(er,RHOstarM);    
    dpdrhor = (prsp-prsm)/(2.00*HTTOL*Wa_n.rho);
    dpdrho  = 0.50*(dpdrhol + dpdrhor);
  }
  cAvg = sqrt(Wa_n.p*dpde/(Wa_n.rho*Wa_n.rho)+dpdrho + 2.0*Wa_n.k/3.0 + (2.0*Wa_n.k*dpde)/(3.0*Wa_n.rho));   

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl_n.lambda_x();
  lambdas_r = Wr_n.lambda_x();
  //using Glaister's computed sound-speed
  lambdas_a = Wa_n.lambda_x(cAvg);
  //lambdas_a = Wa_n.lambda_x(); - with ROE !!

  //--------------------------end Glaister stuff

    /* Determine the intermediate state flux. */

    wavespeed.x = min(lambdas_l[1], lambdas_a[1]);
    wavespeed.y = max(lambdas_r[4], lambdas_a[4]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    return (wavespeed);

}

/**********************************************************************
 * Routine: FluxHLLL_x (Timur Linde's flux function, x-direction)     *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * Linde approximation for the fluxes.  See Linde (1998).             *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLL(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
  HighTemp2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux, dFrl, dUrl, dFwave;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();

  // Determine the intermediate state flux.
  wavespeed_l = min(lambdas_l[1],lambdas_a[1]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
  
  if (wavespeed_l >= ZERO) {
    Flux = Wl.Fx();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.Fx();
  } else {
    dUrl = Wr.U() - Wl.U();
    dFrl = Wr.Fx() - Wl.Fx();
    wavespeed_m = Wa.v.x;
    da = Wa.rho;
    ca = Wa.c();
    dU = (fabs(dUrl.rho)/da + 
	  fabs(dUrl.dv.x)/(da*ca) + 
	  fabs(dUrl.dv.y)/(da*ca) + 
	  fabs(dUrl.E)/(da*ca*ca));
    if (dU <= TOLER) {
      alpha = ZERO;
    } else {
      dU = ONE/dU;
      dFwave = dFrl - wavespeed_m*dUrl;
      alpha = ONE - (fabs(dFwave.rho)/(da*ca) + 
		     fabs(dFwave.dv.x)/(da*ca*ca) + 
		     fabs(dFwave.dv.y)/(da*ca*ca) + 
		     fabs(dFwave.E)/(da*ca*ca*ca))*dU;
      alpha = max(ZERO,alpha);
    }
    
    Flux = ((wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
	    +(wavespeed_l*wavespeed_r)*
	    (ONE-(ONE-max(wavespeed_m/wavespeed_r,
			  wavespeed_m/wavespeed_l))*alpha)*dUrl)/
      (wavespeed_r-wavespeed_l);
  }

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxHLLL(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxHLLL(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxHLLL_n (Timur Linde's flux function, n-direction)     *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the Linde approximation to specify the intermediate  *
 * state flux in terms of the rotated solution states.  See Linde     *
 * (1998).                                                            *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLL_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxHLLL(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxHLLL_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxHLLL_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxGHLLL_x (Timur Linde's flux function, x-direction     *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * Linde approximation for the fluxes.  See Linde (1998).             *
 *  *** Modified for HIGHTEMP cases, with Glaister Averages           *
 **********************************************************************/
HighTemp2D_cState FluxGHLLL(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
  HighTemp2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  HighTemp2D_cState Flux, dFrl, dUrl, dFwave;

   //Glaister-algorithm definitions
  int deltaE, deltaRHO;
  double pspl, psml, pspr, psmr, prl, plr, plsp, plsm, prsp, prsm;
  double de, drho, dpde, dpdrho, dpder, dpdel, dpdrhol, dpdrhor;
  double RHOstarP, RHOstarM, EstarP, EstarM;
  double hl, hr, el, er, srhol, srhor, ha, ea, cAvg;
 
  // Evaluate the modified Roe-average (Glaister) primitive solution state.

  // Determine the left and right state specific enthalpies and square
  // roots of the density.
  el    = Wl.e();
  er    = Wr.e();
  hl    = el + Wl.p/Wl.rho;
  hr    = er + Wr.p/Wr.rho;
  srhol = sqrt(Wl.rho);
  srhor = sqrt(Wr.rho);
  
  // Determine the appropriate Roe averages.
  Wa.rho   = srhol*srhor;
  Wa.v.x   = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
  Wa.v.y   = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
  Wa.k     = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
  Wa.omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
  ha = (srhol*hl + srhor*hr)/(srhol+srhor);
  // Glaister modification
  ea = (srhol*el + srhor*er)/(srhol+srhor);
  Wa.p = Wa.rho*(ha - ea);
 
  //calculation of modified sound speed for Average State  
  RHOstarP = Wa.rho*HTONEPLUST;
  RHOstarM = Wa.rho*HTONEMINT;
  EstarP = ea*HTONEPLUST;
  EstarM = ea*HTONEMINT;
  de = (er-el);
  drho = (Wr.rho-Wl.rho);
 
  //Find dpde derivative
  if (fabs(de) >= (TOLER3*ea)) {
    prl = Tgas_p(er, Wl.rho);
    plr = Tgas_p(el, Wr.rho);
    dpde = 0.5*(Wr.p - plr + prl - Wl.p)/de;
  }else{
    pspl  = Tgas_p(EstarP,Wl.rho);
    psml  = Tgas_p(EstarM,Wl.rho);
    dpdel = (pspl-psml)/(2.00*HTTOL*ea);       
    pspr  = Tgas_p(EstarP,Wr.rho);
    psmr  = Tgas_p(EstarM,Wr.rho);
    dpder = (pspr-psmr)/(2.00*HTTOL*ea);    
    dpde  = 0.50*(dpdel+dpder);
  }

  //Find dpdrho derivative
  if (fabs(drho) >= (TOLER3*Wa.rho)) {
    prl = Tgas_p(er, Wl.rho);
    plr = Tgas_p(el, Wr.rho);
    dpdrho = 0.5*(Wr.p - prl + plr - Wl.p)/drho;
  }else{
    plsp    = Tgas_p(el,RHOstarP);
    plsm    = Tgas_p(el,RHOstarM); 
    dpdrhol = (plsp-plsm)/(2.00*HTTOL*Wa.rho);
    prsp    = Tgas_p(er,RHOstarP);
    prsm    = Tgas_p(er,RHOstarM);    
    dpdrhor = (prsp-prsm)/(2.00*HTTOL*Wa.rho);
    dpdrho  = 0.50*(dpdrhol + dpdrhor);
  }
  cAvg = sqrt(Wa.p*dpde/(Wa.rho*Wa.rho)+dpdrho + 2.0*Wa.k/3.0 + (2.0*Wa.k*dpde)/(3.0*Wa.rho));   

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  //using Glaister's computed sound-speed
  lambdas_a = Wa.lambda_x(cAvg);

  // Determine the intermediate state flux.
  wavespeed_l = min(lambdas_l[1],lambdas_a[1]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
  
  if (wavespeed_l >= ZERO) {
    Flux = Wl.Fx();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.Fx();
  } else {
    dUrl = Wr.U() - Wl.U();
    dFrl = Wr.Fx() - Wl.Fx();
    wavespeed_m = Wa.v.x;
    da = Wa.rho;
    //ca = Wa.c();
    ca = cAvg;
    dU = (fabs(dUrl.rho)/da + 
	  fabs(dUrl.dv.x)/(da*ca) + 
	  fabs(dUrl.dv.y)/(da*ca) + 
	  fabs(dUrl.E)/(da*ca*ca));
    if (dU <= TOLER) {
      alpha = ZERO;
    } else {
      dU = ONE/dU;
      dFwave = dFrl - wavespeed_m*dUrl;
      alpha = ONE - (fabs(dFwave.rho)/(da*ca) + 
		     fabs(dFwave.dv.x)/(da*ca*ca) + 
		     fabs(dFwave.dv.y)/(da*ca*ca) + 
		     fabs(dFwave.E)/(da*ca*ca*ca))*dU;
      alpha = max(ZERO,alpha);
    }
    
    Flux = ((wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
	    +(wavespeed_l*wavespeed_r)*
	    (ONE-(ONE-max(wavespeed_m/wavespeed_r,
			  wavespeed_m/wavespeed_l))*alpha)*dUrl)/
      (wavespeed_r-wavespeed_l);
  }

  // Return the solution flux.
  return Flux;
}

HighTemp2D_cState FluxGHLLL(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxGHLLL(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxGHLLL_n (Timur Linde's flux function, n-direction)    *
 *        *** modified for HighTemp with Glaister's Avgs              *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the Linde approximation to specify the intermediate  *
 * state flux in terms of the rotated solution states.  See Linde     *
 * (1998).                                                            *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxGHLLL_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxGHLLL(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxGHLLL_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxGHLLL_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxHLLC (Tito Toro's flux function, x-direction)         *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the HLLC *
 * approximation for the fluxes.  See Toro, Spruce, and Speares       *
 * (Shock Waves Vol.4 1994).                                          *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLC(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m;
  double al, ar, CL, CR, Z, ql, qr;
  double um, pm, aml, amr;
  
  HighTemp2D_cState Flux, Uml, Umr;

  // Determine the left, intermediate, and right wave speeds.
  al  = Wl.a();
  ar  = Wr.a();
  CL  = Wl.v.x+TWO*al/Wl.gm1;
  CR  = Wr.v.x-TWO*ar/Wr.gm1;
  Z   = (ar/al)*pow((Wl.p/Wr.p),HALF*Wl.gm1/Wl.g);
  um  = (CL*Z+CR)/(ONE+Z);
  aml = al-HALF*Wl.gm1*(um-Wl.v.x);
  pm  = Wl.p*pow((aml/al),TWO*Wl.g/Wl.gm1);
  amr = ar+HALF*Wr.gm1*(um-Wr.v.x);
  pm  = HALF*(pm + Wr.p*pow((amr/ar),TWO*Wr.g/Wr.gm1));

  if (pm/Wl.p <= ONE) {
    ql = ONE;
  } else {
    ql = sqrt(ONE+HALF*((Wl.g+ONE)/Wl.g)*(pm/Wl.p));
  }
  wavespeed_l = Wl.v.x - ql*al;
  
  if (pm/Wr.p <= ONE) {
    qr = ONE;
  } else {
    qr = sqrt(ONE+HALF*((Wr.g+ONE)/Wr.g)*(pm/Wr.p));
  }
  wavespeed_r = Wr.v.x + qr*ar;
  
  wavespeed_m = um;
  
  // Determine the intermediate state flux.
  if (wavespeed_l >= ZERO) {
    Flux = Wl.Fx();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.Fx();
  } else  if (wavespeed_m >= ZERO) {
    Uml = HighTemp2D_cState(Wl.rho,
			 Wl.rho*wavespeed_m,
			 Wl.rho*Wl.v.y,
			 Wl.E()+Wl.rho*(wavespeed_m-Wl.v.x)*(wavespeed_m+Wl.p/(Wl.rho*(wavespeed_l-Wl.v.x)))
			 )*((wavespeed_l-Wl.v.x)/(wavespeed_l-wavespeed_m));
    Flux = Wl.Fx()+wavespeed_l*(Uml-Wl.U());
  } else {
    Umr = HighTemp2D_cState(Wr.rho,
			 Wr.rho*wavespeed_m,
			 Wr.rho*Wr.v.y,
			 Wr.E()+Wr.rho*(wavespeed_m-Wr.v.x)*(wavespeed_m+Wr.p/(Wr.rho*(wavespeed_r-Wr.v.x)))
			 )*((wavespeed_r-Wr.v.x)/(wavespeed_r-wavespeed_m));
    Flux = Wr.Fx()+wavespeed_r*(Umr-Wr.U());
  }

  // Return solution flux.
  return Flux;

}

HighTemp2D_cState FluxHLLC(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxHLLC(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxHLLC_n (Tito Toro's flux function, n-direction)       *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the HLLC approximation to specify the intermediate   *
 * state in terms of the rotated solution states.  See Toro, Spruce,  *
 * and Speares (Shock Waves Vol.4 1994).                              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxHLLC_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxHLLC(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

HighTemp2D_cState FluxHLLC_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxHLLC_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxVanLeer (Bram Van Leer's flux vector splitting flux   *
 *                       function, x-direction)                       *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the flux *
 * vector splitting approximation of Bram Van Leer.                   *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxVanLeer(const HighTemp2D_pState &Wl,
				  const HighTemp2D_pState &Wr) {

  HighTemp2D_cState F, Fl, Fr;
  double M, Mp, Mn;

  M = Wl.v.x/Wl.a();
  if (M < -ONE) {
    Fl.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mp =  0.25*Wl.rho*Wl.a()*sqr(ONE + M);
    Fl[1] = Mp;
    Fl[2] = Mp*((TWO*Wl.a()/Wl.g)*(HALF*Wl.gm1*M + ONE));
    Fl[3] = Mp*Wl.v.y;
    Fl[4] = Mp*((TWO*Wl.a2()/(Wl.g*Wl.g-ONE))*sqr(HALF*Wl.gm1*M+ONE) + HALF*sqr(Wl.v.y));
  } else if (M > ONE) {
    Fl[1] = Wl.rho*Wl.v.x;
    Fl[2] = Wl.rho*sqr(Wl.v.x) + Wl.p;
    Fl[3] = Wl.rho*Wl.v.x*Wl.v.y;
    Fl[4] = Wl.v.x*Wl.H();
  }

  M = Wr.v.x/Wr.a();
  if (-M < -ONE) {
    Fr.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mn = -0.25*Wr.rho*Wr.a()*sqr(ONE - M);
    Fr[1] = Mn;
    Fr[2] = Mn*((TWO*Wr.a()/Wr.g)*(HALF*Wr.gm1*M - ONE));
    Fr[3] = Mn*Wr.v.y;
    Fr[4] = Mn*((TWO*Wr.a2()/(Wr.g*Wl.g-ONE))*sqr(HALF*Wr.gm1*M-ONE) + HALF*sqr(Wr.v.y));
  } else if (-M > ONE) {
    Fr[1] = Wr.rho*Wr.v.x;
    Fr[2] = Wr.rho*sqr(Wr.v.x) + Wr.p;
    Fr[3] = Wr.rho*Wr.v.x*Wr.v.y;
    Fr[4] = Wr.v.x*Wr.H();
  }

  // Return solution flux.
  return Fl + Fr;

}

HighTemp2D_cState FluxVanLeer(const HighTemp2D_cState &Ul,
				  const HighTemp2D_cState &Ur) {
  return FluxVanLeer(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxVanLeer (Bram Van Leer's flux vector splitting flux   *
 *                       function, n-direction)                       *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the flux *
 * vector splitting approximation of Bram Van Leer.                   *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxVanLeer_n(const HighTemp2D_pState &Wl,
				    const HighTemp2D_pState &Wr,
				    const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated, Vl_rotated, Vr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxVanLeer(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

HighTemp2D_cState FluxVanLeer_n(const HighTemp2D_cState &Ul,
				    const HighTemp2D_cState &Ur,
				    const Vector2D &norm_dir) {
  return FluxVanLeer_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxVanLeer_MB (Bram Van Leer's flux vector splitting     *
 *                          flux function, x-direction)               *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the flux *
 * vector splitting approximation of Bram Van Leer for a moving frame *
 * of reference (see Chassaing, Gerolymos, and Vallet. AIAA Journal   *
 * Vol. 41 No. 10 2003).                                              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxVanLeer_MB(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr,
				     const Vector2D &V) {

  HighTemp2D_cState Fl, Fr;
  double M, Mp, Mn;

  M = (Wl.v.x - V.x)/Wl.a();
  if (M < -ONE) {
    Fl.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mp = 0.25*(M + ONE)*(M + ONE);
    Fl[1] = Mp*Wl.rho*Wl.a();
    Fl[2] = Mp*(Wl.rho*Wl.a()*Wl.v.x + (TWO - M)*Wl.p);
    Fl[3] = Mp*Wl.rho*Wl.a()*Wl.v.y;
    Fl[4] = Mp*(Wl.a()*Wl.E() + (TWO - M)*Wl.v.x*Wl.p);
  } else if (M > ONE) {
    Fl[1] = Wl.rho*(Wl.v.x - V.x);
    Fl[2] = Wl.rho*Wl.v.x*(Wl.v.x - V.x) + Wl.p;
    Fl[3] = Wl.rho*Wl.v.y*(Wl.v.x - V.x);
    Fl[4] = Wl.E()*(Wl.v.x - V.x) + Wl.v.x*Wl.p;
  }

  M = (Wr.v.x - V.x)/Wr.a();
  if (-M < -ONE) {
    Fr.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mn = -0.25*(M - ONE)*(M - ONE);
    Fr[1] = Mn*Wr.rho*Wr.a();
    Fr[2] = Mn*(Wr.rho*Wr.a()*Wr.v.x - (TWO + M)*Wr.p);
    Fr[3] = Mn*Wr.rho*Wr.a()*Wr.v.y;
    Fr[4] = Mn*(Wr.a()*Wr.E() - (TWO + M)*Wr.v.x*Wr.p);
  } else if (-M > ONE) {
    Fr[1] = Wr.rho*(Wr.v.x - V.x);
    Fr[2] = Wr.rho*Wr.v.x*(Wr.v.x - V.x) + Wr.p;
    Fr[3] = Wr.rho*Wr.v.y*(Wr.v.x - V.x);
    Fr[4] = Wr.E()*(Wr.v.x - V.x) + Wr.v.x*Wr.p;
  }

  // Return solution flux.
  return Fl + Fr;

}

/**********************************************************************
 * Routine: FluxVanLeer_MB (Bram Van Leer's flux vector splitting     *
 *                          flux function, n-direction)               *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the flux *
 * vector splitting approximation of Bram Van Leer for a moving frame *
 * of reference (see Chassaing, Gerolymos, and Vallet. AIAA Journal   *
 * Vol. 41 No. 10 2003).                                              *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxVanLeer_MB_n(const HighTemp2D_pState &Wl,
				       const HighTemp2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_pState Wl_nonstationary, Wr_nonstationary;
  HighTemp2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Move left and right states back into nonstationary frame.
  Wl_nonstationary = Wl; Wl_nonstationary.v += V;
  Wr_nonstationary = Wr; Wr_nonstationary.v += V;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl_nonstationary,norm_dir);
  Wr_rotated.Rotate(Wr_nonstationary,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxVanLeer_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: FluxAUSM (M.-S. Liou and C.J. Steffen Jr.'s Advection     *
 *                    Upstream Splitting Method flux function,        *
 *                    x-direction)                                    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the AUSM *
 * (Advection Upstream Splitting Method) approximation for the        *
 * fluxes.  See M.-S. Liou and C. J. Steffen, Jr (J. Comp. Physics    *
 * Vol. 107 1993).                                                    *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxAUSM(const HighTemp2D_pState &Wl,
			       const HighTemp2D_pState &Wr) {

  HighTemp2D_cState Flux;
  double al, ar, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;

  // Determine the left and right state sound speed and Mach numbers:
  al = Wl.a(); Ml = Wl.v.x/al;
  ar = Wr.a(); Mr = Wr.v.x/ar;

  // Determine the left state split Mach number:
  if (fabs(Ml) <= ONE) {
    Mplus = 0.25*sqr(Ml+ONE);
    pplus = 0.25*Wl.p*sqr(Ml+ONE)*(TWO-Ml);
  } else {
    Mplus = 0.5*(Ml + fabs(Ml));
    pplus = 0.5*Wl.p*(Ml + fabs(Ml))/Ml;
  }

  // Determine the right state split Mach number:
  if (fabs(Mr) <= ONE) {
    Mminus = -0.25*sqr(Mr-ONE);
    pminus = 0.25*Wr.p*sqr(Mr-ONE)*(TWO+Mr);
  } else {
    Mminus = 0.5*(Mr - fabs(Mr));
    pminus = 0.5*Wr.p*(Mr - fabs(Mr))/Mr;
  }

  // Determine the intermediate state Mach number and pressure:
  Mhalf = Mplus + Mminus;
  phalf = pplus + pminus;

  // Determine the intermediate state solution convective flux:
  Flux[1] = 0.5*Mhalf*(Wl.rho*al + Wr.rho*ar);
  Flux[2] = 0.5*Mhalf*(Wl.rho*al*Wl.v.x + Wr.rho*ar*Wr.v.x);
  Flux[3] = 0.5*Mhalf*(Wl.rho*al*Wl.v.y + Wr.rho*ar*Wr.v.y);
  Flux[4] = 0.5*Mhalf*(al*Wl.H() + ar*Wr.H());

  // Add the numerical dissipation to the intermediate state solution flux:
  Flux[1] -= 0.5*fabs(Mhalf)*(Wr.rho*ar - Wl.rho*al);
  Flux[2] -= 0.5*fabs(Mhalf)*(Wr.rho*ar*Wr.v.x - Wl.rho*al*Wl.v.x);
  Flux[3] -= 0.5*fabs(Mhalf)*(Wr.rho*ar*Wr.v.y - Wl.rho*al*Wl.v.y);
  Flux[4] -= 0.5*fabs(Mhalf)*(ar*Wr.H() - al*Wl.H());

  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;

  // Return solution flux.
  return Flux;

}

HighTemp2D_cState FluxAUSM(const HighTemp2D_cState &Ul,
			       const HighTemp2D_cState &Ur) {
  return FluxAUSM(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxAUSM_n (M.-S. Liou and C.J. Steffen Jr.'s Advection   *
 *                      Upstream Splitting Method flux function,      *
 *                      n-direction)                                  *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the AUSM approximation to specify the intermediate   *
 * state in terms of the rotated solution states.  See M.-S. Liou and *
 * C. J. Steffen, Jr (J. Comp. Physics Vol. 107 1993).                *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxAUSM_n(const HighTemp2D_pState &Wl,
				 const HighTemp2D_pState &Wr,
				 const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSM(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

HighTemp2D_cState FluxAUSM_n(const HighTemp2D_cState &Ul,
				 const HighTemp2D_cState &Ur,
				 const Vector2D &norm_dir) {
  return FluxAUSM_n(Ul.W(),Ur.W(),norm_dir);
}


/**********************************************************************
 * Routine: FluxAUSMplus (Liou's updated Advection Upstream Splitting *
 *                        Method flux function, x-direction)          *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * AUSM+ (updated AUSM scheme) approximation for the fluxes.  See     *
 * M.-S. Liou (NASA-TM-106524, 1994).                                 *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxAUSMplus(const HighTemp2D_pState &Wl,
				   const HighTemp2D_pState &Wr) {

  HighTemp2D_cState Flux;
  double beta = 0.125, alpha = 0.1875;
  double ahalf, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;

  // Determine the intermediate state sound speed:
  ahalf = HALF*(Wl.a() + Wr.a());
  //ahalf = sqrt(Wl.a()*Wr.a());
  //ahalf = min(Wl.a(),Wr.a());

  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;

  // Determine the left state split Mach number:
  if (fabs(Ml) < ONE) {
    Mplus = 0.25*sqr(Ml+ONE) + beta*sqr(Ml*Ml-ONE);
    pplus = (0.25*sqr(Ml+ONE)*(TWO-Ml) + alpha*Ml*sqr(Ml*Ml-ONE))*Wl.p;
  } else {
    Mplus = 0.5*(Ml + fabs(Ml));
    pplus = 0.5*Wl.p*(ONE+double(sgn(Ml)));
  }

  // Determine the right state split Mach number:
  if (fabs(Mr) < ONE) {
    Mminus = -0.25*sqr(Mr-ONE) - beta*sqr(Mr*Mr-ONE);
    pminus = (0.25*sqr(Mr-ONE)*(TWO+Mr) - alpha*Mr*sqr(Mr*Mr-ONE))*Wr.p;
  } else {
    Mminus = 0.5*(Mr - fabs(Mr));
    pminus = 0.5*Wr.p*(ONE-double(sgn(Mr)));
  }

  // Determine the intermediate state Mach number and pressure:
  Mhalf = Mplus + Mminus;
  phalf = pplus + pminus;

  // Determine the intermediate state solution convective flux:
  Flux[1] = 0.5*Mhalf*(Wl.rho + Wr.rho);
  Flux[2] = 0.5*Mhalf*(Wl.rho*Wl.v.x + Wr.rho*Wr.v.x);
  Flux[3] = 0.5*Mhalf*(Wl.rho*Wl.v.y + Wr.rho*Wr.v.y);
  Flux[4] = 0.5*Mhalf*(Wl.H() + Wr.H());

  // Add the numerical dissipation to the intermediate state solution flux:
  Flux[1] -= 0.5*fabs(Mhalf)*(Wr.rho - Wl.rho);
  Flux[2] -= 0.5*fabs(Mhalf)*(Wr.rho*Wr.v.x - Wl.rho*Wl.v.x);
  Flux[3] -= 0.5*fabs(Mhalf)*(Wr.rho*Wr.v.y - Wl.rho*Wl.v.y);
  Flux[4] -= 0.5*fabs(Mhalf)*(Wr.H() - Wl.H());

  // Multiply the intermediate state solution flux by the intermediate
  // state sound speed:
  Flux *= ahalf;

  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;

  // Return solution flux.
  return Flux;

}

HighTemp2D_cState FluxAUSMplus(const HighTemp2D_cState &Ul,
				   const HighTemp2D_cState &Ur) {
  return FluxAUSMplus(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxAUSMplus_n (M.-S. Liou's updated Advection Upstream   *
 *                          Splitting Method flux function,           *
 *                          n-direction)                              *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the AUSM+ approximation to specify the intermediate  *
 * state in terms of the rotated solution states.  See M.-S. Liou     *
 * (NASA-TM-106524, 1994).                                            *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxAUSMplus_n(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr,
				     const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSMplus(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

HighTemp2D_cState FluxAUSMplus_n(const HighTemp2D_cState &Ul,
				     const HighTemp2D_cState &Ur,
				     const Vector2D &norm_dir) {
  return FluxAUSMplus_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxAUSMplus_up (Liou's updated Advection Upstream        * 
 *                           Splitting Method flux function for all   *
 *                           speeds,  x-direction)                    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * AUSM+-up (updated AUSM scheme) approximation for the fluxes.  See  *
 * M.-S. Liou (J. Comp. Physics 2006).                                *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_pState &Wl,
			      const HighTemp2D_pState &Wr) {

  HighTemp2D_cState Flux; //Convected_Quantities;
  double beta = 0.125, sigma = 1.0, Kp =0.25, Ku = 0.5/*0.75*/;
  double alpha, rhohalf, mass_flux_half;
  double ahalf, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;
  //double al, ar, atilde_l, atilde_r;

  ahalf = HALF*(Wl.a() + Wr.a());
  rhohalf = HALF*(Wl.rho + Wr.rho); 
  
  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;

  // Determine the reference Mach number, scaling function and coefficient
  double M2_bar, M2_ref, fa;
  M2_bar = (Wl.v.x*Wl.v.x + Wr.v.x*Wr.v.x)/(TWO*ahalf*ahalf);
  M2_ref = min(ONE, max(M2_bar, Wl.Mref*Wl.Mref));
  if (M2_ref > ONE || M2_ref < 0.0) cout << "\nM2_ref out of range";
  //fa = sqrt(M2_ref)*(TWO - sqrt(M2_ref));
  fa = sqrt(sqr(ONE - M2_ref)*M2_bar + FOUR*M2_ref)/(ONE + M2_ref);
  if (fa > ONE || fa <= ZERO) cout << "\nfa out of range";
  alpha = (3.0/16.0)*(-4.0 + 5.0*fa*fa);
  if (alpha < (-3.0/4.0)  ||  alpha > (3.0/16.0)) cout << "\nalpha out of range";

  // Determine the left state split Mach number:
  if (fabs(Ml) >= ONE) {
    Mplus = Mplus_1(Ml);
    pplus = Mplus_1(Ml)/Ml;
  } else {
    Mplus = Mplus_2(Ml) * (1.0 - 16.0*beta*Mminus_2(Ml));
    pplus = Mplus_2(Ml) * ((2.0 - Ml) - 16.0*alpha*Ml*Mminus_2(Ml));
  }

  // Determine the right state split Mach number:
  if (fabs(Mr) >= ONE) {
    Mminus = Mminus_1(Mr);
    pminus = Mminus_1(Mr)/Mr;        
  } else {
    Mminus = Mminus_2(Mr) * (1.0 + 16.0*beta*Mplus_2(Mr));
    pminus = Mminus_2(Mr) * ((-2.0 - Mr) + 16.0*alpha*Mr*Mplus_2(Mr));
  } 

  // Determine the intermediate state Mach number, pressure and mass flux:
  Mhalf = Mplus + Mminus
    - (Kp/fa)*max((ONE - sigma*M2_bar), ZERO)*(Wr.p - Wl.p)/(rhohalf*ahalf*ahalf);

  phalf = pplus*Wl.p + pminus*Wr.p
    - Ku*pplus*pminus*TWO*rhohalf*(fa*ahalf)*(Wr.v.x - Wl.v.x);

  mass_flux_half = (Mhalf > ZERO) ? ahalf*Mhalf*Wl.rho : ahalf*Mhalf*Wr.rho; 


  // Determine the intermediate state convective solution flux:
  if (mass_flux_half  > ZERO) {
    Flux[1] = ONE;
    Flux[2] = Wl.v.x; 
    Flux[3] = Wl.v.y;
    Flux[4] =  Wl.H()/Wl.rho;
  } else {
    Flux[1] = ONE;
    Flux[2] = Wr.v.x; 
    Flux[3] = Wr.v.y;
    Flux[4] =  Wr.H()/Wr.rho;
  } //end if

  Flux *= mass_flux_half;

  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;

  // Return solution flux.
  return Flux;

}

HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_cState &Ul,
			      const HighTemp2D_cState &Ur) {
  return FluxAUSMplusUP(Ul.W(),Ur.W());
}


/**********************************************************************
 * Routine: FluxAUSMplusUP (Liou's updated Advection Upstream Splitting *
 *                        Method flux function, x-direction)          *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * AUSM+ (updated AUSM scheme) approximation for the fluxes.  See     *
 * M.-S. Liou (Journal of Computational Physics, 2005).               *
 *                                                                    *
 **********************************************************************/
/*
HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_pState &Wl,
				   const HighTemp2D_pState &Wr) {

  HighTemp2D_cState Flux;
  double beta = 0.125, Kp = 0.25, sigma = 1.0, Ku = 0.75;
  double alpha, ahalf, rhohalf, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;
  double Mbar_sq, Mmax, Mo_sq, fa, MAXval, Mp, Pu, mflux, Mref;

  // Determine the intermediate state sound speed:
  ahalf = HALF*(Wl.a() + Wr.a());
  rhohalf = HALF*(Wl.rho + Wr.rho);
  //ahalf = sqrt(Wl.a()*Wr.a());

  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;

  ///to be changed:
  Mref = 0.2;

  Mbar_sq = (sqr(Wl.v.x)+sqr(Wr.v.x))/(2.0*ahalf*ahalf);
  Mmax = max(Mbar_sq,Mref*Mref);
  Mo_sq = min(ONE,Mmax);
  fa = 2.0*sqrt(Mo_sq) - Mo_sq;
  alpha = 0.1875*(-4.0+5.0*sqr(fa));

  // Determine the left state split Mach number:
  if (fabs(Ml) < ONE) {
    Mplus = 0.25*sqr(Ml+ONE) + beta*sqr(Ml*Ml-ONE);
    pplus = (0.25*sqr(Ml+ONE)*(TWO-Ml) + alpha*Ml*sqr(Ml*Ml-ONE))*Wl.p;
  } else {
    Mplus = 0.5*(Ml + fabs(Ml));
    pplus = 0.5*Wl.p*(ONE+fabs(Ml)/Ml);
  }

  // Determine the right state split Mach number:
  if (fabs(Mr) < ONE) {
    Mminus = -0.25*sqr(Mr-ONE) - beta*sqr(Mr*Mr-ONE);
    pminus = (0.25*sqr(Mr-ONE)*(TWO+Mr) - alpha*Mr*sqr(Mr*Mr-ONE))*Wr.p;
  } else {
    Mminus = 0.5*(Mr - fabs(Mr));
    pminus = 0.5*Wr.p*(ONE-fabs(Mr)/Mr);
  }

  //Additional pressure diffusion term to enhance calcs of low Mach flow
  MAXval = max(ONE - sigma*Mbar_sq, 0.0);
  Mp = -(Kp/fa)*MAXval*(Wr.p - Wl.p)/(rhohalf*ahalf*ahalf); 
  
  //additional diffusion term Pu
  Pu = -Ku*(pplus/Wl.p)*(pminus/Wr.p)*(Wl.rho+Wr.rho)*fa*ahalf*(Wr.v.x-Wl.v.x);
  
  // Determine the intermediate state Mach number and pressure:
  Mhalf = Mplus + Mminus + Mp;
  phalf = pplus + pminus + Pu;

  //Define the mass flux
  if (Mhalf>0.0)
    mflux = ahalf*Mhalf*Wl.rho;
  else
    mflux = ahalf*Mhalf*Wr.rho;

  // Determine the intermediate state solution convective flux:
  Flux[1] = mflux;
  if (mflux>0.0) {
    Flux[2] = mflux*Wl.v.x;
    Flux[3] = mflux*Wl.v.y;
    Flux[4] = mflux*Wl.H();
  }else {
    Flux[2] = mflux*Wr.v.x;
    Flux[3] = mflux*Wr.v.y;
    Flux[4] = mflux*Wr.H();
  }

  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;

  // Return solution flux.
  return Flux;
}

HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_cState &Ul,
				   const HighTemp2D_cState &Ur) {
  return FluxAUSMplusUP(Ul.W(),Ur.W());
}
*/


/**********************************************************************
 * Routine: FluxAUSMplus_up_n (M.-S. Liou's Advection Upstream        *
 *                             Splitting Method flux function for     *
 *                             all speeds, n-direction)               *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the AUSM+-up approximation to specify the            * 
 * intermediate  state in terms of the rotated solution states.       *
 * See M.-S. Liou (J. Comp. Physics 2006).                            *
 **********************************************************************/
HighTemp2D_cState FluxAUSMplusUP_n(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr,
				     const Vector2D &norm_dir) {

  HighTemp2D_pState Wl_rotated, Wr_rotated;
  HighTemp2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSMplusUP(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

HighTemp2D_cState FluxAUSMplusUP_n(const HighTemp2D_cState &Ul,
				     const HighTemp2D_cState &Ur,
				     const Vector2D &norm_dir) {
  return FluxAUSMplusUP_n(Ul.W(),Ur.W(),norm_dir);
}


/**********************************************************************
 * Routine: ViscousFlux_n                                             *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
HighTemp2D_cState ViscousFlux_n(const Vector2D &X,
				    HighTemp2D_pState &W,
				    const HighTemp2D_pState &dWdx,
				    const HighTemp2D_pState &dWdy,
				    const Vector2D &norm_dir,
				    int Axisymmetric)
{

  HighTemp2D_cState Gx, Gy, U;

  // Compute the intermediate state viscous stress tensor and heat flux
  // vector.
  W.ComputeViscousTerms(dWdx,dWdy,X,Axisymmetric);
  U = W.U(); U.tau = W.tau; U.q = W.q;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  Gx = U.Gx(dWdx);
  Gy = U.Gy(dWdy);

  // Return the intermediate state solution viscous flux.
  return Gx*norm_dir.x + Gy*norm_dir.y;

}

/**********************************************************************
 * Routine: WallShearStress                                           *
 *                                                                    *
 * This routine returns the wall shear stress in the coordinate frame *
 * defined by the inward unit normal of the wall.  The gradient of    *
 * the velocity field is determined using a Green-Gauss integration   *
 * over the half-diamond-path.                                        *
 *                                                                    *
 **********************************************************************/
double WallShearStress(const HighTemp2D_pState &W1,
		       const Vector2D &X1,
		       const Vector2D &X2,
		       const Vector2D &X3,
		       const Vector2D &nhat) {

  double A;
  Vector2D n21, n32, n13, that;
  HighTemp2D_pState W2, W3, W_face, dWdx, dWdy, dWdt, dWdn;

  // Determine the tangential direction at the wall.
  that = Vector2D(nhat.y,-nhat.x);

  // Initialze W2 and W3.
  W2.Vacuum(); W2.rho = W1.rho; W2.p = W1.p;
  W3.Vacuum(); W3.rho = W1.rho; W3.p = W1.p;

  // Determine the normals of the faces and the area of the region of
  // Green-Gauss integration.
  n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x));
  n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x));
  n13 = Vector2D((X1.y-X3.y),-(X1.x-X3.x));
  A = HALF*((X2-X1)^(X3-X1));

  // Compute Green-Gauss integration on triangle.
  W_face = HALF*(W2+W1);
  dWdx = W_face*n21.x;
  dWdy = W_face*n21.y;
  W_face = HALF*(W3+W2);
  dWdx += W_face*n32.x;
  dWdy += W_face*n32.y;
  W_face = HALF*(W1+W3);
  dWdx += W_face*n13.x;
  dWdy += W_face*n13.y;
  dWdx /= A;
  dWdy /= A;

  // Rotate the gradients of the velocity components into the reference
  // frame of the wall.
  dWdt.v.x = dWdx.v.x*that.x + dWdy.v.x*that.y;
  dWdt.v.y = dWdx.v.y*that.x + dWdy.v.y*that.y;
  dWdn.v.x = dWdx.v.x*nhat.x + dWdy.v.x*nhat.y;
  dWdn.v.y = dWdx.v.y*nhat.x + dWdy.v.y*nhat.y;

  // Return the wall shear stress.
  return (W2.mu()+W2.muT())*(dWdt.v.y + dWdn.v.x);

}

double WallShearStress2(const Vector2D &X,
			const Vector2D &X1,
			const HighTemp2D_pState &W1,
			const HighTemp2D_pState &dW1dx,
			const HighTemp2D_pState &dW1dy,
			const Vector2D &nhat) {

  HighTemp2D_pState W2, dWdx, dWdy, dWds, dWdt, dWdn;
  Vector2D dX, that;
  double ds;

  W2.Vacuum(); W2.rho = W1.rho; W2.p = W1.p;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.

  dX = X-X1; ds = dX.abs(); dX /= ds;

  dWds = (W2-W1)/ds;

  dWdx = dW1dx + (dWds - dW1dx*dX.x)*nhat.x/dot(nhat,dX);
  dWdy = dW1dy + (dWds - dW1dy*dX.y)*nhat.y/dot(nhat,dX);

  // Determine the tangential direction at the wall.
  that = Vector2D(nhat.y,-nhat.x);

  // Rotate the gradients of the velocity components into the reference
  // frame of the wall.
  dWdt.v.x = dWdx.v.x*that.x + dWdy.v.x*that.y;
  dWdt.v.y = dWdx.v.y*that.x + dWdy.v.y*that.y;
  dWdn.v.x = dWdx.v.x*nhat.x + dWdy.v.x*nhat.y;
  dWdn.v.y = dWdx.v.y*nhat.x + dWdy.v.y*nhat.y;

  // Return the wall shear stress.
  return W2.mu()*(dWdt.v.y + dWdn.v.x);

}

/**********************************************************************
 * Routine: ShearStress                                               *
 *                                                                    *
 * This routine returns shear stress in the coordinate frame defined  *
 * by the given normal.                                               *
 *                                                                    *
 **********************************************************************/
double ShearStress(const HighTemp2D_pState &W,
		   const HighTemp2D_pState &dWdx,
		   const HighTemp2D_pState &dWdy,
		   const Vector2D &nhat) {

  Vector2D that;
  HighTemp2D_pState dWdt, dWdn;

  // Determine the tangential direction at the wall.
  that = Vector2D(nhat.y,-nhat.x);

  // Rotate the gradients of the velocity components into the reference
  // frame of the wall.
  dWdt.v.x = dWdx.v.x*that.x + dWdy.v.x*that.y;
  dWdn.v.x = dWdx.v.x*nhat.x + dWdy.v.x*nhat.y;
  dWdt.v.y = dWdx.v.y*that.x + dWdy.v.y*that.y;
  dWdn.v.y = dWdx.v.y*nhat.x + dWdy.v.y*nhat.y;

  // Return the wall shear stress.
  return W.mu()*(dWdt.v.y + dWdn.v.x);

}

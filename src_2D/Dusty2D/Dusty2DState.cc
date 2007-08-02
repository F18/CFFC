/**********************************************************************
 * Dusty2DState.cc                                                    *
 *                                                                    *
 * Subroutines for 2D dusty solution state classes.                   *
 *                                                                    *
 **********************************************************************/

// Include 2D Dusty solution state header file.

#ifndef _DUSTY2D_STATE_INCLUDED
#include "Dusty2DState.h"
#endif // _DUSTY2D_STATE_INCLUDED

/**********************************************************************
 * Dusty2D_pState -- Create storage and assign gas-phase constants.   *
 **********************************************************************/
double Dusty2D_pState::g    = GAMMA_AIR;
double Dusty2D_pState::gm1  = GAMMA_AIR-ONE;
double Dusty2D_pState::gm1i = ONE/(GAMMA_AIR-ONE);
double Dusty2D_pState::R    = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
double Dusty2D_pState::cp   = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double Dusty2D_pState::v1   = AIR_c1;
double Dusty2D_pState::v2   = AIR_c2;
double Dusty2D_pState::v3   = AIR_c3;
double Dusty2D_pState::v4   = AIR_c4;
double Dusty2D_pState::v5   = AIR_c5;
int Dusty2D_pState::flow_type = FLOWTYPE_INVISCID;

/**********************************************************************
 * Dusty2D_cState -- Create storage and assign gas-phase constants.   *
 **********************************************************************/
double Dusty2D_cState::g    = GAMMA_AIR;
double Dusty2D_cState::gm1  = GAMMA_AIR-ONE;
double Dusty2D_cState::gm1i = ONE/(GAMMA_AIR-ONE);
double Dusty2D_cState::R    = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
double Dusty2D_cState::cp   = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double Dusty2D_cState::v1   = AIR_c1;
double Dusty2D_cState::v2   = AIR_c2;
double Dusty2D_cState::v3   = AIR_c3;
double Dusty2D_cState::v4   = AIR_c4;
double Dusty2D_cState::v5   = AIR_c5;
int Dusty2D_cState::flow_type = FLOWTYPE_INVISCID;

/**********************************************************************
 * Dusty2D_pState -- Create storage and assign dust-phase constants.  *
 **********************************************************************/
int Dusty2D_pState::NUM_CMP_PART = 0;
int Dusty2D_pState::NUM_VAR_PART = 0;
int Dusty2D_pState::drag_law = DRAG_LAW_HAIDER_LEVENSPIEL;
double Dusty2D_pState::d1 = 0.1806;
double Dusty2D_pState::d2 = 0.6459;
double Dusty2D_pState::d3 = 0.4251;
double Dusty2D_pState::d4 = 6880.05;
double Dusty2D_pState::d5 = -1.0;
double Dusty2D_pState::cm   = CM_APHTPB;
double Dusty2D_pState::dp   = DP_APHTPB;
double Dusty2D_pState::rhop = RHOP_APHTPB;
double Dusty2D_pState::mp   = MP_APHTPB;
double Dusty2D_pState::qe   = QE_APHTPB;

/**********************************************************************
 * Dusty2D_cState -- Create storage and assign dust-phase constants.  *
 **********************************************************************/
int Dusty2D_cState::NUM_CMP_PART = 0;
int Dusty2D_cState::NUM_VAR_PART = 0;
int Dusty2D_cState::drag_law = DRAG_LAW_HAIDER_LEVENSPIEL;
double Dusty2D_cState::d1 = 0.1806;
double Dusty2D_cState::d2 = 0.6459;
double Dusty2D_cState::d3 = 0.4251;
double Dusty2D_cState::d4 = 6880.05;
double Dusty2D_cState::d5 = -1.0;
double Dusty2D_cState::cm   = CM_APHTPB;
double Dusty2D_cState::dp   = DP_APHTPB;
double Dusty2D_cState::rhop = RHOP_APHTPB;
double Dusty2D_cState::mp   = MP_APHTPB;
double Dusty2D_cState::qe   = QE_APHTPB;

/**********************************************************************
 * Dusty2D_pState -- Create storage and assign turbulent static       *
 *                   variables.                                       *
 **********************************************************************/
// Turbulent boundary-layer constants:
double Dusty2D_pState::yplus_o = 10.0;
double Dusty2D_pState::C = 5.0;
double Dusty2D_pState::von_karman = 0.41;
double Dusty2D_pState::yplus_sublayer = 5.0;
double Dusty2D_pState::yplus_buffer_layer = 30.0;
double Dusty2D_pState::yplus_outer_layer = 100.0;
// k-omega closure coefficients:
double Dusty2D_pState::PrT = 0.90;
double Dusty2D_pState::Cmu = 0.090;
double Dusty2D_pState::beta_k_o = 0.09;
double Dusty2D_pState::beta_omega_o = 0.072;
double Dusty2D_pState::sigma_k = 0.50;
double Dusty2D_pState::sigma_omega = 0.50;
double Dusty2D_pState::alpha = 0.52;
double Dusty2D_pState::xi = 1.50;
double Dusty2D_pState::Mto = 0.25;

/**********************************************************************
 * Dusty2D_cState -- Create storage and assign turbulent static       *
 *                   variables.                                       *
 **********************************************************************/
// Turbulent boundary-layer constants:
double Dusty2D_cState::yplus_o = 10.0;
double Dusty2D_cState::C = 5.0;
double Dusty2D_cState::von_karman = 0.41;
double Dusty2D_cState::yplus_sublayer = 5.0;
double Dusty2D_cState::yplus_buffer_layer = 30.0;
double Dusty2D_cState::yplus_outer_layer = 100.0;
// k-omega coefficients:
double Dusty2D_cState::PrT = 0.90;
double Dusty2D_cState::Cmu = 0.090;
double Dusty2D_cState::beta_k_o = 0.09;
double Dusty2D_cState::beta_omega_o = 0.072;
double Dusty2D_cState::sigma_k = 0.50;
double Dusty2D_cState::sigma_omega = 0.50;
double Dusty2D_cState::alpha = 0.52;
double Dusty2D_cState::xi = 1.50;
double Dusty2D_cState::Mto = 0.25;

/**********************************************************************
 * Dusty2D_pState -- Create storage and assign propellant constants.  *
 **********************************************************************/
double Dusty2D_pState::rhos   = RHOS_APHTPB;
double Dusty2D_pState::n      = N_APHTPB;
double Dusty2D_pState::beta   = BETA_APHTPB;
double Dusty2D_pState::Tf     = TF_APHTPB;
double Dusty2D_pState::Ts     = TS_APHTPB;
double Dusty2D_pState::alphas = ALPHAS_APHTPB;

/**********************************************************************
 * Dusty2D_cState -- Create storage and assign propellant constants.  *
 **********************************************************************/
double Dusty2D_cState::rhos   = RHOS_APHTPB;
double Dusty2D_cState::n      = N_APHTPB;
double Dusty2D_cState::beta   = BETA_APHTPB;
double Dusty2D_cState::Tf     = TF_APHTPB;
double Dusty2D_cState::Ts     = TS_APHTPB;
double Dusty2D_cState::alphas = ALPHAS_APHTPB;

/**********************************************************************
 * Dusty2D_pState -- Create storage for the total number of variables.*
 **********************************************************************/
int Dusty2D_pState::NUM_VAR_DUSTY2D = NUM_VAR_BASE;
int* Dusty2D_pState::cindex = NULL;
int* Dusty2D_pState::pindex = NULL;

/**********************************************************************
 * Dusty2D_cState -- Create storage for the total number of variables.*
 **********************************************************************/
int Dusty2D_cState::NUM_VAR_DUSTY2D = NUM_VAR_BASE;
int* Dusty2D_cState::cindex = NULL;
int* Dusty2D_cState::pindex = NULL;

/**********************************************************************
 * Routine: Riemann (Exact Riemann solver, x-direction)               *
 *                                                                    *
 * This function uses a Newton-Raphson interative procedure to obtain *
 * the exact solution to the Riemann problem for the 2D Dusty         *
 * equations in the x-direction, returning the intermediate state     *
 * variables along the ray x/t=0.  See Gottlieb and Groth (1987).     *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState Riemann(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  int number_of_iterations;
  
  double al, ar, CL, CR, Z;
  double dml, dmr, vm, vml, vmr, pm, aml, amr;
  double msl, pml, dpmldum, msr, pmr, dpmrdum;
  double vsl, vhl, vtl, vsr, vhr, vtr;

  Dusty2D_pState W; W.Vacuum();

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
	W = Dusty2D_pState(dml,vm,vml,pm);
	return W;
      } else {
	vm = (Wl.gm1*Wl.v.x + TWO*al)/(Wl.g + ONE);
	pm = Wl.p*pow((vm/al),TWO*Wl.g/Wl.gm1);
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(vm);
	W = Dusty2D_pState(dml,vm,vml,pm);
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
	W = Dusty2D_pState(dmr,vm,vmr,pm);
 	return W;
      } else {
	vm = (Wr.gm1*Wr.v.x - TWO*ar)/(Wr.g + ONE);
	pm = Wr.p*pow((-vm/ar),TWO*Wr.g/Wr.gm1);
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(vm);
	W = Dusty2D_pState(dmr,vm,vmr,pm);
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
  if (vm >= ZERO) {
    if (vm < Wl.v.x) {
      aml = al*sqrt(((Wl.g+ONE)+Wl.gm1*pm/Wl.p)/
	            ((Wl.g+ONE)+Wl.gm1*Wl.p/pm));
      vsl = Wl.v.x+msl*al;
      if (vsl >= ZERO) {
	return Dusty2D_pState(Wl);                
      } else {
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(aml);
	return Dusty2D_pState(dml,vm,vml,pm);
      }
    } else {
      vhl = Wl.v.x-al;
      vtl = vm-aml;
      if (vhl >= ZERO) {
	return Wl;
      } else if (vtl <= ZERO) {
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(aml);
	W = Dusty2D_pState(dml,vm,vml,pm);
	return W;
      } else {
	vm = (Wl.gm1*Wl.v.x+TWO*al)/(Wl.g+ONE);
	pm = Wl.p*pow((vm/al),TWO*Wl.g/Wl.gm1);
	vml = Wl.v.y;
	dml = Wl.g*pm/sqr(vm);
	W = Dusty2D_pState(dml,vm,vml,pm);
	return W;
      }
    }
  } else {
    if (vm > Wr.v.x) {
      amr = ar*sqrt(((Wr.g+ONE)+Wr.gm1*pm/Wr.p)/
		    ((Wr.g+ONE)+Wr.gm1*Wr.p/pm));
      vsr = Wr.v.x+msr*ar;
      if (vsr <= ZERO) {
	return Wr;
      } else {
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(amr);
	W = Dusty2D_pState(dmr,vm,vmr,pm);
	return W;
      }
    } else {
      vhr = Wr.v.x + ar;
      vtr = vm + amr;
      if (vhr <= ZERO) {
	return Dusty2D_pState(Wr);
      } else if (vtr >= ZERO) {
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(amr);
	W = Dusty2D_pState(dmr,vm,vmr,pm);
	return W;
      } else {
	vm = (Wr.gm1*Wr.v.x-TWO*ar)/(Wr.g+ONE);
	pm = Wr.p*pow((-vm/ar),TWO*Wr.g/Wr.gm1);
	vmr = Wr.v.y;
	dmr = Wr.g*pm/sqr(vm);
	W = Dusty2D_pState(dmr,vm,vmr,pm);
	return W;
      }
    }
  }

}

/**********************************************************************
 * Routine: RoeAverage (Roe Averages)                                 *
 *                                                                    *
 * This function returns the Roe-averaged primitive solution state    *
 * given left and right primitive solution variables. See Roe (1981). *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState RoeAverage(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {
  
  double hl, hr, srhol, srhor, aa2, ha;
  Dusty2D_pState Wa; Wa.rho = ZERO; Wa.p = ZERO;//Wa.Vacuum();

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
  aa2 = Wa.gm1*(ha - HALF*(sqr(Wa.v.x) + sqr(Wa.v.y)) - Wa.k);
  Wa.p = Wa.rho*aa2/Wa.g;

  // Return the Roe-averged state.
  return Wa;

}

/**********************************************************************
 * Routine: Rotate                                                    *
 *                                                                    *
 * This function returns the solution in the local rotated frame.     *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState Rotate(const Dusty2D_pState &W, const Vector2D &norm_dir) {

  Dusty2D_pState W_rotated(W);
  double cos_angle = norm_dir.x;
  double sin_angle = norm_dir.y;
  W_rotated.rho   = W.rho;
  W_rotated.v.x   =   W.v.x*cos_angle + W.v.y*sin_angle;
  W_rotated.v.y   = - W.v.x*sin_angle + W.v.y*cos_angle;
  W_rotated.p     = W.p;
  W_rotated.k     = W.k;
  W_rotated.omega = W.omega;
  if (W.NUM_VAR_PART) W_rotated.Wp = Rotate(W.Wp,norm_dir);
  // Return the rotated state.
  return W_rotated;

}

Dusty2D_cState Rotate(const Dusty2D_cState &U, const Vector2D &norm_dir) {

  Dusty2D_cState U_rotated(U);
  double cos_angle = norm_dir.x;
  double sin_angle = norm_dir.y;
  U_rotated.rho    = U.rho;
  U_rotated.dv.x   =   U.dv.x*cos_angle + U.dv.y*sin_angle;
  U_rotated.dv.y   = - U.dv.x*sin_angle + U.dv.y*cos_angle;
  U_rotated.E      = U.E;
  U_rotated.dk     = U.dk;
  U_rotated.domega = U.domega;
  if (U.NUM_VAR_PART) U_rotated.Up = Rotate(U.Up,norm_dir);
  // Return the rotated state.
  return U_rotated;

}

/**********************************************************************
 * Routine: RotationMatrix                                            *
 *                                                                    *
 * This function returns the rotation matrix.                         *
 *                                                                    *
 **********************************************************************/
extern DenseMatrix RotationMatrix(const Dusty2D_pState &Wdum,
				  const Vector2D &norm_dir) {

  double cos_angle = norm_dir.x; 
  double sin_angle = norm_dir.y;
    
  DenseMatrix Omega(Wdum.NUM_VAR_DUSTY2D,Wdum.NUM_VAR_DUSTY2D);
  Omega.identity();
  // Gas-phase variables.
  Omega(1,1) = cos_angle;
  Omega(1,2) = sin_angle;
  Omega(2,1) = -sin_angle;
  Omega(2,2) = cos_angle;
  // Particle-phase variables.
  // Return the rotation matrix.
  return Omega;

}

/**********************************************************************
 * Routine: iRotationMatrix                                           *
 *                                                                    *
 * This function returns the inverse of the rotation matrix.          *
 *                                                                    *
 **********************************************************************/
extern DenseMatrix iRotationMatrix(const Dusty2D_pState &Wdum,
				   const Vector2D &norm_dir) {

  double cos_angle = norm_dir.x; 
  double sin_angle = norm_dir.y;
    
  DenseMatrix Omega(Wdum.NUM_VAR_DUSTY2D,Wdum.NUM_VAR_DUSTY2D);
  Omega.identity();
  // Gas-phase variables.
  Omega(1,1) =  cos_angle;
  Omega(1,2) = -sin_angle;
  Omega(2,1) =  sin_angle;
  Omega(2,2) =  cos_angle;
  // Particle-phase variables.
  // Return the inverse of the rotation matrix.
  return Omega;

}

/**********************************************************************
 * Routine: Translate                                                 *
 *                                                                    *
 * This function returns the solution in a stationary frame.          *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState Translate(const Dusty2D_pState &W, const Vector2D &V) {

  Dusty2D_pState W_translated;

  W_translated.rho = W.rho;
  W_translated.v.x = W.v.x - V.x;
  W_translated.v.y = W.v.y - V.y;
  W_translated.p   = W.p;
  W_translated.k   = W.k;
  W_translated.omega   = W.omega;
  if (W.NUM_VAR_PART) W_translated.Wp = Translate(W.Wp,V);
  // Return the translated state.
  return W_translated;

}

/**********************************************************************
 * Routine: Reflect                                                   *
 *                                                                    *
 * This function returns the reflected solution state in a given      *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
// Dusty2D_pState Reflect(const Dusty2D_pState &W,
// 		       const Vector2D &norm_dir,
// 		       const int &reflect_flag) {

//   Dusty2D_pState Wr, Wn;

//   // Apply the frame rotation and calculate the primitive solution 
//   // state variables in the local rotated frame defined by the unit 
//   // normal vector.
//   Wr = Rotate(W,norm_dir);

//   // Reflect the gas-phase normal velocity in the rotated frame.
//   Wr.v.x = -Wr.v.x;

//   // Rotate back to the original Cartesian reference frame.
//   Wn = Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

//   // The action of the reflection boundary condition is dependent on the
//   // formulation type.  If the single-velocity component formulation is
//   // used then the reflection boundary condition is calculated normally.
//   // If the multi-velocity formulation is being used then the
//   // reflection boundary condition acts as normal reflection for a flux
//   // calculation or as a constant extrapolation boundary condition for
//   // filling the ghost cells.
// //   if (W.NUM_VAR_PART) {
// //     if (W.NUM_CMP_PART == PARTICLE2D_SINGLE_VELOCITY_FORMULATION ||
// //  	reflect_flag) Wn.Wp = Reflect(W.Wp,norm_dir);
// //     else Wn.Wp = W.Wp;
// //   }
//   if (W.NUM_VAR_PART) {
//       Wn.Wp = Reflect(W.Wp,norm_dir);
// //     if (W.NUM_CMP_PART == PARTICLE2D_SINGLE_VELOCITY_FORMULATION ||
// //  	reflect_flag) {
// //       Wn.Wp = Reflect(W.Wp,norm_dir);
// //     } else {
// //       Wn.Wp = W.Wp;
// // //       for (int nc = 0; nc < W.Wp.NUM_CMP_PARTICLE2D; nc++) {
// // // 	if (W.Wp[nc].sigma < NANO) {
// // // 	} else {
// // // 	  if (dot(norm_dir,W.Wp[nc].u) > ZERO) Wn.Wp[nc] = W.Wp[nc];
// // // 	}
// // //       }
// //     }
//   }

//   // Return the reflected state.
//   return Wn;

// }
Dusty2D_pState Reflect(const Dusty2D_pState &W,
		       const Vector2D &norm_dir) {

  Dusty2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr = Rotate(W,norm_dir);

  // Reflect the gas-phase normal velocity in the rotated frame.
  Wr.v.x = - Wr.v.x;

  // Rotate back to the original Cartesian reference frame.
  Wn = Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

  // The action of the reflection boundary condition is dependent on the
  // formulation type.  If the single-velocity component formulation is
  // used then the reflection boundary condition is calculated normally.
  // If the multi-velocity formulation is being used then the reflection
  // boundary condition acts as normal reflection for a flux calculation
  // or as a constant extrapolation boundary condition for filling the
  // ghost cells.
  if (W.NUM_VAR_PART) Wn.Wp = Reflect(W.Wp,norm_dir);

  // Return the reflected state.
  return Wn;

}

/**********************************************************************
 * Routine: Absorb                                                    *
 *                                                                    *
 * This function returns the absorbed (constant extrapolation for     *
 * dusty phase and reflection for the gas phase) solution state in a  * 
 * given direction given the primitive solution variables and the     *
 * unit normal vector in the direction of interest.                   *
 *                                                                    *
 **********************************************************************/
// Dusty2D_pState Absorb(const Dusty2D_pState &W,
// 		      const Vector2D &norm_dir,
// 		      const int &reflect_flag) {
  
//   Dusty2D_pState Wn;

//   // Apply the reflection boundary condition.
//   Wn = Reflect(W,norm_dir,reflect_flag);

//   // Absorb the particle-phase.
//   if (W.NUM_VAR_PART) Wn.Wp = Absorb(W.Wp,norm_dir);

//   // Return the absorbed state.
//   return Wn;

// }
Dusty2D_pState Absorb(const Dusty2D_pState &W,
		      const Vector2D &norm_dir) {
  
  Dusty2D_pState Wn;

  // Apply the reflection boundary condition.
  Wn = Reflect(W,norm_dir);

  // Absorb the particle-phase.
  if (W.NUM_VAR_PART) Wn.Wp = Absorb(W.Wp,norm_dir);

  // Return the absorbed state.
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
Dusty2D_pState WallViscousHeatFlux(const Dusty2D_pState &W,
				   const Vector2D &norm_dir) {

  return MovingWallHeatFlux(W,norm_dir,ZERO);

}

/**********************************************************************
 * Routine: WallViscousIsothermal                                     *
 *                                                                    *
 * This function returns the no-slip solution state with a fixed wall *
 * temperature in a given direction given the primitive solution      *
 * variables and the unit normal vector in the direction of interest. *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState WallViscousIsothermal(const Dusty2D_pState &W,
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
Dusty2D_pState MovingWallHeatFlux(const Dusty2D_pState &W,
				  const Vector2D &norm_dir,
				  const double &Vwall) {

  Dusty2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr = Rotate(W,norm_dir);

  // Use the Roe-averaged value to find the correct tangential velocity.
  Wr.v.x = - Wr.v.x;
  Wr.v.y = - TWO*Vwall - Wr.v.y;
  if (W.NUM_VAR_PART) {
    for (int np = 0; np < W.NUM_CMP_PART; np++) {
      if (Wr.Wp[np].sigma > ZERO) {
 	Wr.Wp[np].u.x = - Wr.Wp[np].u.x;
 	Wr.Wp[np].u.y = - TWO*Vwall - Wr.Wp[np].u.y;
      }
    }
  }

  // Apply the adiabatic viscous wall boundary condition to the
  // turbulent variables.
  Wr.k = ZERO;

  // Rotate back to the original Cartesian reference frame.
  Wn = Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_pState MovingWallIsothermal(const Dusty2D_pState &W,
				    const Vector2D &norm_dir,
				    const double &Vwall,
				    const double &Twall) {

  Dusty2D_pState Wr, Wn;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr = Rotate(W,norm_dir);

  // Use the Roe-averaged value to find the correct tangential velocity.
  Wr.v.x = - Wr.v.x;
  Wr.v.y = - TWO*Vwall - Wr.v.y;
  if (W.NUM_VAR_PART) {
    for (int np = 0; np < W.NUM_CMP_PART; np++) {
      if (Wr.Wp[np].sigma > ZERO) {
	Wr.Wp[np].u.x = - Wr.Wp[np].u.x;
	Wr.Wp[np].u.y = - TWO*Vwall - Wr.Wp[np].u.y;
      }
    }
  }

  // Apply the specified wall temperature.
  Wr.rho = W.p/(W.R*Twall);
  if (W.NUM_VAR_PART) {
    for (int np = 0; np < W.NUM_CMP_PART; np++) {
      if (Wr.Wp[np].sigma > ZERO) Wr.Wp[np].Tp = Twall;
    }
  }

  // Apply the isothermal viscous wall boundary condition to the
  // turbulent variables.
  Wr.k = ZERO;

  // Rotate back to the original Cartesian reference frame.
  Wn = Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_pState BurningSurface(const Dusty2D_pState &W, const Vector2D &norm_dir) {

  double g = W.g, R = W.R, vxr, vyr, pr, ar, alphas; 
  double sigmabs, rhobs, vxbs, vxbsp, rbs = ZERO, rbsp, vybs, pbs, pbso, ars;
  double cos_angle, sin_angle, vx, vy, C1, C2;
  Dusty2D_pState Wbs;
  int iteration = 0;

  if (!W.NUM_VAR_PART) alphas = ZERO;
  else alphas = W.alphas;

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
    pbs   = pbso - (((ONE-alphas)*W.rhos*rbs - rhobs*vxbs)/
		    ((ONE-alphas)*W.rhos*rbsp - rhobs*(vxbs/pbso + vxbsp)));
    if (pbs <= ZERO) pbs = HALF*pbso;
    iteration++;
  } while (fabs(ONE - (pbso/pbs)) > TOLER/HUNDRED);
  if (fabs(vxbs) < TOLER) {
    sigmabs = ZERO;
  } else {
    sigmabs = alphas*W.rhos*rbs/vxbs;
  }
  vybs = ZERO;
  vxbs = -vxbs; // Required to correct for burning direction.

  // Rotate back to the original reference frame.
  Wbs.rho = rhobs;
  Wbs.v.x = vxbs*cos_angle - vybs*sin_angle;
  Wbs.v.y = vxbs*sin_angle + vybs*cos_angle;
  Wbs.p   = pbs;
  if (Wbs.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    double sigmav = 0.035;
    double lw = 0.000200;
    Wbs.k     = sqr(sigmav)*sqr(vxbs);
    Wbs.omega = pow(W.Cmu,-0.75)/(lw*Wbs.beta_k_o);
  }
  if (W.NUM_VAR_PART) Wbs.Wp = Particle2D_pComponents(Particle2D_pState(sigmabs,Wbs.v.x,Wbs.v.y,W.Tf));

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
Dusty2D_pState RinglebFlow(const Dusty2D_pState &Wdum,
			   const Vector2D &X) {
  double q, k;
  return RinglebFlow(Wdum,X,q,k);
}

Dusty2D_pState RinglebFlow(const Dusty2D_pState &Wdum,
			   const Vector2D &X,
			   double &q, double &k) {

  Dusty2D_pState W; W.Vacuum();
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
  if (k > 5.0/3.0) { W.Standard_Atmosphere(); return W; }//cout << "k = " << k << " > 5/3 @ " << X << endl;
  sin_theta = max(ZERO,min(ONE,q/k));
  theta = TWO*PI-asin(sin_theta);
  sin_theta = sin(theta);
  cos_theta = cos(theta);
  W.rho *= rhoo;
  W.v.x = sqrt(g*po/rhoo)*q*cos_theta;
  if (X.y < ZERO) W.v.x *= -ONE;
  W.v.y = sqrt(g*po/rhoo)*q*sin_theta;
  W.p = po*(W.rho/rhoo)*c*c;

  // Return W state.
  return W;

}

/**********************************************************************
 * Routine: RinglebFlowAverageState                                   *
 *                                                                    *
 * This function returns the average exact solution to Ringleb's flow *
 * for a quadrilateral cell defined by points X1, X2, X3, and X4      *
 * which are defined in a counter-clockwise direction starting in the *
 * south-west corner.  The average solution is determined using a     *
 * ?????th-order numerical quadrature.                                *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState RinglebFlowAverageState(const Dusty2D_pState &Wdum,
				       const Vector2D &Y1,
				       const Vector2D &Y2,
				       const Vector2D &Y3,
				       const Vector2D &Y4) {

  Vector2D X, X1, X2, X3, X4;
  double epsilon1 = -ONE, epsilon2 = ONE, epsilon3 = -ONE, epsilon4 = ONE;
  double eta1 = -ONE, eta2 = -ONE, eta3 = ONE, eta4 = ONE;
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
  Dusty2D_pState W; W.Vacuum();

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
Dusty2D_pState ViscousChannelFlow(const Dusty2D_pState &Wdum,
				  const Vector2D X,
				  const Vector2D Vwall,
				  const double dp,
				  const double length,
				  const double height) {
  Dusty2D_pState W;
  // Compute the exact laminar channel solution.  Note that the
  // pressure, p, and density, rho, must be set before the x-direction
  // velocity component, v.x, since the v.x depends on the dynamic
  // viscosity, mu, which depends on the temperature, T.
  W.rho = DENSITY_STDATM;
  W.p   = PRESSURE_STDATM - dp*(ONE - X.x/length);
  W.v.x = (HALF/W.mu())*dp*X.y*(X.y - height)/length + Vwall.x*X.y/height;
  W.v.y = ZERO;
  // Return W state.
  return W;
}

Dusty2D_pState ViscousChannelFlowVelocity(const Dusty2D_pState &Wdum,
					  const Vector2D X,
					  const Vector2D Vwall,
					  const double dp,
					  const double length,
					  const double height) {
  Dusty2D_pState W;
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
 * Routine: ViscousChannelFlowDerivatives                             *
 *                                                                    *
 * This function will return the exact laminar channel flow solution  *
 * (Couette or Poiseuille flows) given an (x,y)-coordinate where x =  *
 * [0,length] and y = [0,height], an upper wall speed Vwall, and the  *
 * imposed pressure change, dp.  Use ViscousPipeFlow for axisymmetric *
 * problems.  The exact laminar solution can be used as an initial    *
 * solution for a turbulent channel flow.                             *
 *                                                                    *
 **********************************************************************/
void ViscousChannelFlowDerivatives(const Dusty2D_pState &W,
				   Dusty2D_pState &dWdx,
				   Dusty2D_pState &dWdy,
				   const Vector2D X,
				   const Vector2D Vwall,
				   const double dp,
				   const double length,
				   const double height) {

  dWdx.rho = ZERO;
  dWdx.v.x = ZERO;
  dWdx.v.y = ZERO;
  dWdx.p   = dp/length;

  dWdy.rho = ZERO;
  dWdy.v.x = (HALF/W.mu())*dp*(TWO*X.y - height)/length + Vwall.x/height;
  dWdy.v.y = ZERO;
  dWdy.p   = ZERO;

}

/**********************************************************************
 * Routine: ViscousPipeFlow                                           *
 *                                                                    *
 * This function will return the exact laminar flow solution for an   *
 * axisymmetric pipe flow, given the imposed pressure drop, dp, and   *
 * an (r,z)-coordinate where r = [0,radius] and z = [0,length].       *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState ViscousPipeFlow(const Dusty2D_pState &Wdum,
			       const Vector2D X,
			       const double dp,
			       const double length,
			       const double radius) {

  Dusty2D_pState W;
  // Compute the exact laminar pipe solution.  Note that the pressure,
  // p, and density, rho, must be set before the x-direction velocity 
  // component, v.x, since the v.x depends on the dynamic viscosity, 
  // mu, which depends on the temperature, T.
  W.rho = DENSITY_STDATM;
  W.p   = PRESSURE_STDATM - (dp*length)*(ONE - X.x/length);
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
Dusty2D_pState TurbulentPipeFlow(const Dusty2D_pState &Wo,
				 const Vector2D X,
				 const double dp,
				 const double length,
				 const double radius) {

  Dusty2D_pState W;
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
Dusty2D_pState FlatPlate(const Dusty2D_pState &Winf,
			 const Vector2D &X,
			 const double &plate_length,
			 double &eta,
			 double &f,
			 double &fp,
			 double &fpp) {

  Dusty2D_pState W;
  double fo, dn, k1, k2, k3, k4;

  // Initialize variables.
  W.rho = Winf.rho; W.p = Winf.p;
  eta = ZERO; f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
  //dn = 0.000005;
  dn = 0.0005;

  // Return upstream conditions before flat plate, including the leading edge.
  if (X.x <= ZERO || X.y > HALF*plate_length) return Winf;
  //if (X.x < NANO || X.x > plate_length || X.y > plate_length) return Winf;

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(Winf.v.x/(X.x*Winf.nu()));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    W.v.x = fp*Winf.v.x;
    W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/max(X.x,TOLER))*(eta*fp-f);
    return W;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  for (double n = ZERO; n < eta; n += dn) {

    // Store the solution at the start of the iteration.
    fo = f;

    // Increment f:
    k1 = dn*fp;
    k2 = dn*(fp + k1/2.0);
    k3 = dn*(fp + k2/2.0);
    k4 = dn*(fp + k3);
    f = f + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fp:
    k1 = dn*fpp;
    k2 = dn*(fpp + k1/2.0);
    k3 = dn*(fpp + k2/2.0);
    k4 = dn*(fpp + k3);
    fp = fp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fpp:
    k1 = -dn*fo*fpp/2.0;
    k2 = -dn*(fo + dn/2.0)*(fpp + k1/2.0)/2.0;
    k3 = -dn*(fo + dn/2.0)*(fpp + k2/2.0)/2.0;
    k4 = -dn*(fo + dn)*(fpp + k3)/2.0;
    fpp = fpp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

  }

  // Compute the velocity vector at point X.
  W.v.x = fp*Winf.v.x;
  W.v.y = HALF*sqrt(Winf.nu()*Winf.v.x/max(X.x,TOLER))*(eta*fp-f);

  // Return the Blasius solution state.
  return W;

}

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
Dusty2D_pState DrivenCavityFlow(const Dusty2D_pState &Wo,
				const double &l,
				const double &Re) {
  assert(fabs(Re-100.0) < TOLER || fabs(Re-400.0) < TOLER);
  assert(fabs(l-0.001) < TOLER);
  Dusty2D_pState W;
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
Dusty2D_pState BackwardFacingStep(const Dusty2D_pState &Wo,
				  const Vector2D &X,
				  const double &h,
				  const double &ho,
				  const double &Re,
				  const double &M) {
  assert(fabs(Re-100.0) < TOLER || fabs(Re-389.0) < TOLER || fabs(Re-1000.0) < TOLER);
  assert(fabs(M-0.20) < TOLER);
  Dusty2D_pState W;
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
Dusty2D_pState BC_Characteristic(const Dusty2D_pState &Wi,
                                 const Dusty2D_pState &Wo,
	      	                 const Vector2D &norm_dir) {

  Dusty2D_pState Wi_rotated, Wo_rotated, We, We_rotated;
  char pattern;
  double mi, poi, pai, pbi, mab, mac1, mac2, mbc1, mc1c2, mbd, mad;
  double de, ue, ve, pe, ae, ue_rotated, ve_rotated;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by
  // the unit normal vector.
  Wi_rotated = Rotate(Wi,norm_dir);
  Wo_rotated = Rotate(Wo,norm_dir);

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
    We_rotated.rho = Wi_rotated.rho;
    We_rotated.v.x = Wi_rotated.v.x;
    We_rotated.v.y = Wi_rotated.v.y;
    We_rotated.p   = Wi_rotated.p;
    break;        
  case 'b' :
    We_rotated.rho = Wi_rotated.rho*pow(poi,ONE/Wi_rotated.g);
    We_rotated.v.x = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y = Wi_rotated.v.y;
    We_rotated.p   = Wo_rotated.p;
    break;
  case 'c' :
    We_rotated.rho = Wo_rotated.rho;
    We_rotated.v.x = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(poi,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y = Wo_rotated.v.y;
    We_rotated.p   = Wo_rotated.p;
    break;
  case 'd' :
    We_rotated.p   = Wi_rotated.p*pow((TWO/(Wi_rotated.g+ONE)+(Wi_rotated.gm1/(Wi_rotated.g+ONE))*mi),TWO*Wi_rotated.g*Wi_rotated.gm1i);
    We_rotated.v.x = Wi_rotated.v.x - TWO*Wi_rotated.a()*Wi_rotated.gm1i*(pow(We_rotated.p/Wi_rotated.p,HALF*Wi_rotated.gm1/Wi_rotated.g)-ONE);
    We_rotated.v.y = Wi_rotated.v.y;
    We_rotated.rho = Wi_rotated.rho*pow(We_rotated.p/Wi_rotated.p,ONE/Wi_rotated.g);
    break;
  default:
    We_rotated.rho = Wo_rotated.rho;
    We_rotated.v.x = Wo_rotated.v.x;
    We_rotated.v.y = Wo_rotated.v.y;
    We_rotated.p   = Wo_rotated.p;
    break;
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  We = Rotate(We_rotated,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_pState BC_Characteristic_Pressure(const Dusty2D_pState &Wi,
                                          const Dusty2D_pState &Wo,
	      	                          const Vector2D &norm_dir) {

  Dusty2D_pState Wi_rotated, Wo_rotated, Wb, Wb_rotated;
  double mi, ab;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by 
  // the unit normal vector.
  Wi_rotated = Rotate(Wi,norm_dir);
  Wo_rotated = Rotate(Wo,norm_dir);

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
    Wb_rotated.p     = Wo_rotated.p;
    Wb_rotated.rho   = Wi_rotated.rho*pow(Wb_rotated.p/Wi_rotated.p,ONE/Wi_rotated.g);
    ab = sqrt(Wi_rotated.g*Wb_rotated.p/Wb_rotated.rho);
    Wb_rotated.v.x   = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wi_rotated.gm1i;
    Wb_rotated.v.y   = Wi_rotated.v.y;
    Wb_rotated.k     = Wi_rotated.k;
    Wb_rotated.omega = Wi_rotated.omega;

  // Boundary condition for subsonic inflow.  Pressure specified.
  } else if (mi >= -ONE) {
    Wb_rotated.p     = Wo_rotated.p;
    Wb_rotated.rho   = Wo_rotated.rho;
    ab = sqrt(Wo_rotated.g*Wb_rotated.p/Wb_rotated.rho);
    Wb_rotated.v.x   = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*Wo_rotated.gm1i;
    Wb_rotated.v.y   = Wo_rotated.v.y;
    Wb_rotated.k     = Wo_rotated.k;
    Wb_rotated.omega = Wo_rotated.omega;

  // Boundary condition for supersonic inflow.
  } else {
    Wb_rotated.rho   = Wo_rotated.rho;
    Wb_rotated.v.x   = Wo_rotated.v.x;
    Wb_rotated.v.y   = Wo_rotated.v.y;
    Wb_rotated.p     = Wo_rotated.p;
    Wb_rotated.k     = Wo_rotated.k;
    Wb_rotated.omega = Wo_rotated.omega;
  }

  // If particle-phase x-velocity is greater than zero then use the
  // internal state or if less than zero then use the external state.
  if (Wi.NUM_VAR_PART) {
    for (int np = 0; np < Wi.NUM_CMP_PART; np++) {
      if (Wi_rotated.Wp[np].sigma > ZERO) {
	if (Wi_rotated.Wp[np].u.x >= ZERO) {
	  Wb_rotated.Wp[np] = Wi_rotated.Wp[np];
	} else {
	  Wb_rotated.Wp[np] = Wo_rotated.Wp[np];
	}
	if (Wb_rotated.Wp[np].sigma < TOLER || Wb_rotated.Wp[np].Tp < TOLER) {
	  Wb_rotated.Wp[np].Vacuum();
	}
      }
    }
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  Wb = Rotate(Wb_rotated,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_pState BC_Characteristic_Mach_Number(const Dusty2D_pState &Wi,
                                             const Dusty2D_pState &Wo,
	      	                             const Vector2D &norm_dir) {

  Dusty2D_pState Wi_rotated, Wo_rotated, Wb, Wb_rotated;
  double mi, mb, ab;

  // Apply the frame rotation and evaluate interior and imposed
  // boundary solution states in the local rotated frame defined by 
  // the unit normal vector.
  Wi_rotated = Rotate(Wi,norm_dir);
  Wo_rotated = Rotate(Wo,norm_dir);

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

  // If particle-phase x-velocity is greater than zero then use the
  // internal state or if less than zero then use the external state.
  if (Wi.NUM_VAR_PART) {
    for (int np = 0; np < Wi.NUM_CMP_PART; np++) {
      if (Wi_rotated.Wp[np].sigma > ZERO) {
	if (Wi_rotated.Wp[np].u.x >= ZERO) {
	  Wb_rotated.Wp[np] = Wi_rotated.Wp[np];
	} else {
	  Wb_rotated.Wp[np] = Wo_rotated.Wp[np];
	}
	if (Wb_rotated.Wp[np].sigma < TOLER || Wb_rotated.Wp[np].Tp < TOLER) {
	  Wb_rotated.Wp[np].Vacuum();
	}
      }
    }
  }

  // Rotate the resulting boundary state back to the original 
  // Cartesian reference frame.
  Wb = Rotate(Wb_rotated,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_pState WaveSpeedPos(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wsp;
  for (int n = 0; n < NUM_VAR_BASE; n++)
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
Dusty2D_pState WaveSpeedNeg(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wsn;
  for (int n = 0; n < NUM_VAR_BASE; n++)
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
Dusty2D_pState WaveSpeedAbs(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wsa;
  for (int n = 0; n < NUM_VAR_BASE; n++)
    Wsa[n] = fabs(lambdas_a[n]);
  return Wsa;
}

/**********************************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the positive parts of the corrected          *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState HartenFixPos(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wfp;
  Wfp[1] = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  Wfp[2] = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
  Wfp[3] = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
  Wfp[4] = HartenFixPos(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  Wfp[5] = HALF*(lambdas_a[5]+fabs(lambdas_a[5]));
  Wfp[6] = HALF*(lambdas_a[6]+fabs(lambdas_a[6]));
  return Wfp;
}

/**********************************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the negative parts of the corrected          *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState HartenFixNeg(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wfn;
  Wfn[1] = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  Wfn[2] = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
  Wfn[3] = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
  Wfn[4] = HartenFixNeg(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  Wfn[5] = HALF*(lambdas_a[5]-fabs(lambdas_a[5]));
  Wfn[6] = HALF*(lambdas_a[6]-fabs(lambdas_a[6]));
  return Wfn;
}

/**********************************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)                         *
 *                                                                    *
 * This function returns the absolute values of the corrected         *
 * elemental wave speeds or eigenvalues according to the entropy fix  *
 * of Harten (1983).                                                  *
 *                                                                    *
 **********************************************************************/
Dusty2D_pState HartenFixAbs(const Dusty2D_pState &lambdas_a,
                            const Dusty2D_pState &lambdas_l,
                            const Dusty2D_pState &lambdas_r) {
  Dusty2D_pState Wfa;
  Wfa[1] = HartenFixAbs(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  Wfa[2] = fabs(lambdas_a[2]);
  Wfa[3] = fabs(lambdas_a[3]);
  Wfa[4] = HartenFixAbs(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  Wfa[5] = fabs(lambdas_a[5]);
  Wfa[6] = fabs(lambdas_a[6]);
  return Wfa;
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
Dusty2D_cState FluxGodunov_n(const Dusty2D_pState &Wl,
	      	             const Dusty2D_pState &Wr,
                             const Vector2D &norm_dir) {

  Dusty2D_pState W, Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;
  
  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);
  
  // Solve the Riemann problem in the rotated frame and evaluate the 
  // intermediate solution state and flux.
  W = Riemann(Wl_rotated,Wr_rotated);
  Flux_rotated = W.F();
  
  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

Dusty2D_cState FluxGodunov_n(const Dusty2D_cState &Ul,
			     const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxGodunov_MB_n(const Dusty2D_pState &Wl,
				const Dusty2D_pState &Wr,
				const Vector2D &V,
				const Vector2D &norm_dir) {

  Dusty2D_pState W_rotated, Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Solve the Riemann problem in the rotated frame and evaluate the 
  // intermediate state solution flux.
  W_rotated = Riemann(Wl_rotated,Wr_rotated);

  // Transform back into the moving body frame of reference.
  W_rotated.v.x += V_rotated.x;

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = W_rotated.F(V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

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
Dusty2D_cState FluxRoe(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  Dusty2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_BASE : NUM_VAR_BASE-2;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);
//   cout << endl << " Wa =" << Wa;

  // Evaluate the jumps in the primitive solution states.
  dWrl = Wr - Wl;
//   cout << endl << " dWrl =" << dWrl;

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();
//   cout << endl << " lambdas_l =" << lambdas_l;
//   cout << endl << " lambdas_r =" << lambdas_r;
//   cout << endl << " lambdas_a =" << lambdas_a;
//   cout << endl << " HartenFixNeg =" << HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
//   cout << endl << " HartenFixPos =" << HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
//   for (int i = 1; i <= NumVar; i++) {
//     cout << endl << i << endl << Wa.lp_x(i) << endl << Wa.rc_x(i);
//   }

  // Determine the intermediate state flux.
  if (Wa.v.x >= ZERO) {
    Flux = Wl.F();
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] < ZERO) {
	Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
      }
    }
  } else {
    Flux = Wr.F();
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

Dusty2D_cState FluxRoe(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxRoe_n(const Dusty2D_pState &Wl,
	      	         const Dusty2D_pState &Wr,
                         const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRoe(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Error checking: clear particle terms.
  if (Flux.NUM_VAR_PART) Flux.Up.Vacuum();

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxRoe_n(const Dusty2D_cState &Ul,
			 const Dusty2D_cState &Ur,
			 const Vector2D &norm_dir) {
  return FluxRoe_n(Ul.W(),Ur.W(),norm_dir);
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
Dusty2D_cState FluxRoe_MB(const Dusty2D_pState &Wl,
			  const Dusty2D_pState &Wr,
			  const Vector2D &V) {

  Dusty2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_BASE : NUM_VAR_BASE-2;

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
    Flux = Wl.F(V);
    wavespeeds = HartenFixNeg(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] < ZERO)
	Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  } else {
    Flux = Wr.F(V);
    wavespeeds = HartenFixPos(lambdas_a,lambdas_l,lambdas_r);
    for (int i = 1; i <= NumVar; i++) {
      if (wavespeeds[i] > ZERO)
	Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
    }
  }

  // Return the solution flux.
  return Flux;
    
}

Dusty2D_cState FluxRoe_MB(const Dusty2D_cState &Ul,
			  const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxRoe_MB_n(const Dusty2D_pState &Wl,
			    const Dusty2D_pState &Wr,
			    const Vector2D &V,
			    const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_pState Wl_nonstationary, Wr_nonstationary;
  Dusty2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Move left and right states back into nonstationary frame.
  Wl_nonstationary = Translate(Wl,-V);
  Wr_nonstationary = Translate(Wr,-V);

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl_nonstationary,norm_dir);
  Wr_rotated = Rotate(Wr_nonstationary,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRoe_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxRoe_MB_n(const Dusty2D_cState &Ul,
			    const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxRusanov(const Dusty2D_pState &Wl,
			   const Dusty2D_pState &Wr) {

  double wavespeed_max;
  Dusty2D_pState Wa, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux, dUrl;
  int NumVar = (Wl.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA)
               ? NUM_VAR_BASE : NUM_VAR_BASE-2;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U()-Wl.U(); if (dUrl.NUM_CMP_PART) dUrl.Up.Vacuum();

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x();
  lambdas_r = Wr.lambda_x();
  lambdas_a = Wa.lambda_x();

  // Determine the intermediate state flux.
  Flux = HALF*(Wl.F()+Wr.F());
  wavespeeds = HartenFixAbs(lambdas_a,lambdas_l,lambdas_r);

  wavespeed_max = wavespeeds[1];
  for (int i = 2; i <= NumVar; i++)
    wavespeed_max = max(wavespeed_max,wavespeeds[i]);
  
  Flux -= HALF*wavespeed_max*dUrl;

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxRusanov(const Dusty2D_cState &Ul,
			   const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxRusanov_n(const Dusty2D_pState &Wl,
	      	             const Dusty2D_pState &Wr,
                             const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxRusanov(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxRusanov_n(const Dusty2D_cState &Ul,
			     const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxHLLE(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  double wavespeed_l, wavespeed_r;
  Dusty2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux, dUrl;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U() - Wl.U(); if (dUrl.NUM_CMP_PART) dUrl.Up.Vacuum();

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
    Flux = Wl.F();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F();
  } else {
    Flux = (((wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
	     + (wavespeed_l*wavespeed_r)*dUrl)/
	    (wavespeed_r-wavespeed_l));
  }

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLE(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxHLLE_n(const Dusty2D_pState &Wl,
	      	          const Dusty2D_pState &Wr,
                          const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;
  
  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated
  // frame.
  Flux_rotated = FluxHLLE(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and
  // return the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLE_n(const Dusty2D_cState &Ul,
			  const Dusty2D_cState &Ur,
			  const Vector2D &norm_dir) {
  return FluxHLLE_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxHLLE_MB (HLLE's flux function)                        *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxHLLE_MB(const Dusty2D_pState &Wl,
			   const Dusty2D_pState &Wr,
			   const Vector2D &V) {

  double wavespeed_l, wavespeed_r;
  Dusty2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux, dUrl;

  // Evaluate the Roe-average primitive solution state.
  Wa = RoeAverage(Wl,Wr);

  // Evaluate the jumps in the conserved solution states.
  dUrl = Wr.U() - Wl.U();

  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl.lambda_x(V);
  lambdas_r = Wr.lambda_x(V);
  lambdas_a = Wa.lambda_x(V);

  // Determine the intermediate state flux.
  wavespeed_l = min(lambdas_l[1],lambdas_a[1]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
 
  if (wavespeed_l >= ZERO) {
    Flux = Wl.F(V);
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F(V);
  } else {
    Flux = ((wavespeed_r*Wl.F(V)-wavespeed_l*Wr.F(V)) +
	    (wavespeed_l*wavespeed_r)*dUrl)/(wavespeed_r-wavespeed_l);
  }

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLE_MB(const Dusty2D_cState &Ul,
			   const Dusty2D_cState &Ur,
			   const Vector2D &V) {
  return FluxHLLE_MB(Ul.W(),Ur.W(),V);
}

/**********************************************************************
 * Routine: FluxHLLE_MB_n (HLLE's flux function, n-direction)         *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxHLLE_MB_n(const Dusty2D_pState &Wl,
			     const Dusty2D_pState &Wr,
			     const Vector2D &V,
			     const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_pState Wl_nonstationary, Wr_nonstationary;
  Dusty2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Move left and right states back into nonstationary frame.
  Wl_nonstationary = Translate(Wl,-V);
  Wr_nonstationary = Translate(Wr,-V);

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl_nonstationary,norm_dir);
  Wr_rotated = Rotate(Wr_nonstationary,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxHLLE_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLE_MB_n(const Dusty2D_cState &Ul,
			     const Dusty2D_cState &Ur,
			     const Vector2D &V,
			     const Vector2D &norm_dir) {
  return FluxHLLE_MB_n(Ul.W(),Ur.W(),V,norm_dir);
}

/**********************************************************************
 * Routine: HLLE_wavespeeds                                           *
 *                                                                    *
 * This function returns the two acoustic wavespeeds for an arbitrary *
 * direction defined by a unit normal vector in the direction of      *
 * interest, given left and right primitive solution states.          *
 *                                                                    *
 * Note: wavespeed.x = wavespeed_l = lambda minus.                    *
 *       wavespeed.y = wavespeed_r = lambda plus.                     *
 *                                                                    *
 **********************************************************************/
Vector2D HLLE_wavespeeds(const Dusty2D_pState &Wl,
                         const Dusty2D_pState &Wr,
                         const Vector2D &norm_dir) {

  Vector2D wavespeed;
  Dusty2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  

  // Use rotated values to calculate eignvalues.
  Wl_n = Rotate(Wl,norm_dir);
  Wr_n = Rotate(Wr,norm_dir);

  // Evaluate the Roe-average primitive solution state.
  Wa_n = RoeAverage(Wl_n,Wr_n);
    
  // Evaluate the left, right, and average state eigenvalues.
  lambdas_l = Wl_n.lambda_x();
  lambdas_r = Wr_n.lambda_x();
  lambdas_a = Wa_n.lambda_x();

  // Determine the intermediate state flux.
  wavespeed.x = min(lambdas_l[1],lambdas_a[1]);
  wavespeed.y = max(lambdas_r[4],lambdas_a[4]);

  // Choose the minimum and maximum wavespeeds.
  wavespeed.x = min(wavespeed.x,ZERO);
  wavespeed.y = max(wavespeed.y,ZERO);

  // Return the wavespeeds.
  return wavespeed;

}

/**********************************************************************
 * Routine: FluxHLLL_x (Timur Linde's flux function, x-direction)     *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * Linde approximation for the fluxes.  See Linde (1998).             *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxHLLL(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m, da, ca, dU, alpha;
  Dusty2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
  Dusty2D_cState Flux, dFrl, dUrl, dFwave;

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
    Flux = Wl.F();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F();
  } else {
    dUrl = Wr.U() - Wl.U(); if (dUrl.NUM_CMP_PART) dUrl.Up.Vacuum();
    dFrl = Wr.F() - Wl.F();
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
    
    Flux = ((wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
	    +(wavespeed_l*wavespeed_r)*
	    (ONE-(ONE-max(wavespeed_m/wavespeed_r,
			  wavespeed_m/wavespeed_l))*alpha)*dUrl)/
      (wavespeed_r-wavespeed_l);
  }

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLL(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxHLLL_n(const Dusty2D_pState &Wl,
			  const Dusty2D_pState &Wr,
			  const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxHLLL(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLL_n(const Dusty2D_cState &Ul,
			  const Dusty2D_cState &Ur,
			  const Vector2D &norm_dir) {
  return FluxHLLL_n(Ul.W(),Ur.W(),norm_dir);
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
Dusty2D_cState FluxHLLC(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  double wavespeed_l, wavespeed_r, wavespeed_m;
  double al, ar, CL, CR, Z, ql, qr;
  double um, pm, aml, amr;
  
  Dusty2D_cState Flux, Uml, Umr;

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
    Flux = Wl.F();
  } else if (wavespeed_r <= ZERO) {
    Flux = Wr.F();
  } else  if (wavespeed_m >= ZERO) {
    Uml = Dusty2D_cState(Wl.rho,
			 Wl.rho*wavespeed_m,
			 Wl.rho*Wl.v.y,
			 Wl.E()+Wl.rho*(wavespeed_m-Wl.v.x)*(wavespeed_m+Wl.p/(Wl.rho*(wavespeed_l-Wl.v.x)))
			 )*((wavespeed_l-Wl.v.x)/(wavespeed_l-wavespeed_m));
    Flux = Wl.F()+wavespeed_l*(Uml-Wl.U());
  } else {
    Umr = Dusty2D_cState(Wr.rho,
			 Wr.rho*wavespeed_m,
			 Wr.rho*Wr.v.y,
			 Wr.E()+Wr.rho*(wavespeed_m-Wr.v.x)*(wavespeed_m+Wr.p/(Wr.rho*(wavespeed_r-Wr.v.x)))
			 )*((wavespeed_r-Wr.v.x)/(wavespeed_r-wavespeed_m));
    Flux = Wr.F()+wavespeed_r*(Umr-Wr.U());
  }

  // Return solution flux.
  return Flux;

}

Dusty2D_cState FluxHLLC(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxHLLC_n(const Dusty2D_pState &Wl,
	      	          const Dusty2D_pState &Wr,
                          const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxHLLC(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

Dusty2D_cState FluxHLLC_n(const Dusty2D_cState &Ul,
	      	          const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxVanLeer(const Dusty2D_pState &Wl,
			   const Dusty2D_pState &Wr) {

  Dusty2D_cState F, Fl, Fr;
  double M, Mp, Mn;

  M = Wl.v.x/Wl.a();
  if (M < -ONE) {
    Fl.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mp =  0.25*Wl.rho*Wl.a()*(ONE + M)*(ONE + M);
    Fl[1] = Mp;
    Fl[2] = Mp*((TWO*Wl.a()/Wl.g)*(((Wl.g-ONE)/TWO)*M + ONE));
    Fl[3] = Mp*Wl.v.y;
    Fl[4] = Mp*((TWO*Wl.a2()/(Wl.g*Wl.g-ONE))*(((Wl.g-ONE)/TWO)*M+ONE)*(((Wl.g-ONE)/TWO)*M+ONE) + HALF*Wl.v.y*Wl.v.y);
  } else if (M > ONE) {
    Fl[1] = Wl.rho*Wl.v.x;
    Fl[2] = Wl.rho*Wl.v.x*Wl.v.x + Wl.p;
    Fl[3] = Wl.rho*Wl.v.x*Wl.v.y;
    Fl[4] = Wl.v.x*Wl.H();
  }

  M = Wr.v.x/Wr.a();
  if (-M < -ONE) {
    Fr.Vacuum();
  } else if (fabs(M) <= ONE) {
    Mn = -0.25*Wr.rho*Wr.a()*(ONE - M)*(ONE - M);
    Fr[1] = Mn;
    Fr[2] = Mn*((TWO*Wr.a()/Wr.g)*(((Wr.g-ONE)/TWO)*M - ONE));
    Fr[3] = Mn*Wr.v.y;
    Fr[4] = Mn*((TWO*Wr.a2()/(Wr.g*Wl.g-ONE))*(((Wr.g-ONE)/TWO)*M-ONE)*(((Wr.g-ONE)/TWO)*M-ONE) + HALF*Wr.v.y*Wr.v.y);
  } else if (-M > ONE) {
    Fr[1] = Wr.rho*Wr.v.x;
    Fr[2] = Wr.rho*Wr.v.x*Wr.v.x + Wr.p;
    Fr[3] = Wr.rho*Wr.v.x*Wr.v.y;
    Fr[4] = Wr.v.x*Wr.H();
  }

  // Return solution flux.
  return Fl + Fr;

}

Dusty2D_cState FluxVanLeer(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxVanLeer_n(const Dusty2D_pState &Wl,
			     const Dusty2D_pState &Wr,
			     const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated, Vl_rotated, Vr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxVanLeer(Wl_rotated,Wr_rotated);

  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

Dusty2D_cState FluxVanLeer_n(const Dusty2D_cState &Ul,
			     const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxVanLeer_MB(const Dusty2D_pState &Wl,
			      const Dusty2D_pState &Wr,
			      const Vector2D &V) {

  Dusty2D_cState Fl, Fr;
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
Dusty2D_cState FluxVanLeer_MB_n(const Dusty2D_pState &Wl,
				const Dusty2D_pState &Wr,
				const Vector2D &V,
				const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_pState Wl_nonstationary, Wr_nonstationary;
  Dusty2D_cState Flux, Flux_rotated;
  Vector2D V_rotated;

  // Move left and right states back into nonstationary frame.
  Wl_nonstationary = Wl; Wl_nonstationary.v += V;
  Wr_nonstationary = Wr; Wr_nonstationary.v += V;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl_nonstationary,norm_dir);
  Wr_rotated = Rotate(Wr_nonstationary,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxVanLeer_MB(Wl_rotated,Wr_rotated,V_rotated);

  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

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
Dusty2D_cState FluxAUSM(const Dusty2D_pState &Wl, const Dusty2D_pState &Wr) {

  Dusty2D_cState Flux;
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

Dusty2D_cState FluxAUSM(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxAUSM_n(const Dusty2D_pState &Wl,
	      	          const Dusty2D_pState &Wr,
                          const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSM(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

Dusty2D_cState FluxAUSM_n(const Dusty2D_cState &Ul,
	      	          const Dusty2D_cState &Ur,
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
Dusty2D_cState FluxAUSMplus(const Dusty2D_pState &Wl,
			    const Dusty2D_pState &Wr) {

  Dusty2D_cState Flux;
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

Dusty2D_cState FluxAUSMplus(const Dusty2D_cState &Ul, const Dusty2D_cState &Ur) {
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
Dusty2D_cState FluxAUSMplus_n(const Dusty2D_pState &Wl,
			      const Dusty2D_pState &Wr,
			      const Vector2D &norm_dir) {

  Dusty2D_pState Wl_rotated, Wr_rotated;
  Dusty2D_cState Flux, Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated = Rotate(Wl,norm_dir);
  Wr_rotated = Rotate(Wr,norm_dir);

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSMplus(Wl_rotated,Wr_rotated);
  
  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  Flux = Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.  
  return Flux;

}

Dusty2D_cState FluxAUSMplus_n(const Dusty2D_cState &Ul,
			      const Dusty2D_cState &Ur,
			      const Vector2D &norm_dir) {
  return FluxAUSMplus_n(Ul.W(),Ur.W(),norm_dir);
}

/**********************************************************************
 * Routine: FluxSaurel_n (Saurel's Riemann solver for the dusty-phase *
 *                        of a two-phase flow, n-direction)           *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * particle-phase Riemann solver developed by Saurel, Daniel, and     *
 * Loraud (AIAA J. Vol. 32 No. 6 1992) for the particle-phase of a    *
 * disperse and dilute two-phase flow.                                *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxSaurel_n(const Dusty2D_pState &Wl,
			    const Dusty2D_pState &Wr,
			    const Vector2D &norm_dir) {

  Dusty2D_cState Flux; Flux[1] = ZERO; Flux[4] = ZERO;

  // Evaluate the particle-phase intermediate solution state.
  Flux.Up[0].FluxSaurel_n(Wl.Wp[0],Wr.Wp[0],norm_dir,Wl.cm);

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: FluxSaurel_MB (Saurel's Riemann solver for the dusty      *
 *                         phase of a two-phase flow, n-direction)    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * Riemann solver outlined by Saurel for the dusty phase of a two-    *
 * phase gas flow.  See Saurel, Daniel, and Loraud.  Altered for a    *
 * moving frame of reference.                                         *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxSaurel_MB_n(const Dusty2D_pState &Wl,
			       const Dusty2D_pState &Wr,
			       const Vector2D &V,
			       const Vector2D &norm_dir) {

  Dusty2D_cState Flux; Flux[1] = ZERO; Flux[4] = ZERO;

  // Evaluate the particle-phase intermediate solution state.
  Flux.Up[0].FluxSaurel_MB_n(Wl.Wp[0],Wr.Wp[0],V,norm_dir,Wl.cm);

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: FluxMultiVelocity_n                                       *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by computing the  *
 * flux from the average of the left and right solution states.       *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxMultiVelocity_n(const Dusty2D_pState &Wl,
				   const Dusty2D_pState &Wr,
				   const Vector2D &norm_dir) {

  Dusty2D_cState Flux; Flux[1] = ZERO; Flux[4] = ZERO;

  // Evaluate the particle-phase intermediate solution state.
  Flux.Up = FluxMultiVelocity_n(Wl.Wp,Wr.Wp,norm_dir,Wl.cm);

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: FluxMultiVelocity_MB_n                                    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by computing the  *
 * flux from the average of the left and right solution states.       *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState FluxMultiVelocity_MB_n(const Dusty2D_pState &Wl,
				      const Dusty2D_pState &Wr,
				      const Vector2D &V,
				      const Vector2D &norm_dir) {

  Dusty2D_cState Flux; Flux[1] = ZERO; Flux[4] = ZERO;

  // Evaluate the particle-phase intermediate solution state.
  Flux.Up = FluxMultiVelocity_MB_n(Wl.Wp,Wr.Wp,V,norm_dir,Wl.cm);

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: ViscousFlux_n                                             *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState ViscousFlux_n(const Vector2D &X,
			     Dusty2D_pState &W,
			     const Dusty2D_pState &dWdx,
			     const Dusty2D_pState &dWdy,
			     const Vector2D &norm_dir,
			     const int &Axisymmetric,
			     const int &adiabatic_flag) {

  Dusty2D_cState Gx, Gy, U;

  // Compute the intermediate state viscous stress tensor and heat flux
  // vector.
  W.ComputeViscousTerms(dWdx,dWdy,X,Axisymmetric,adiabatic_flag);
  U = W.U(); U.tau = W.tau; U.q = W.q;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  Gx = U.Gx(dWdx);
  Gy = U.Gy(dWdy);

  // Return the intermediate state solution viscous flux.
  return Gx*norm_dir.x + Gy*norm_dir.y;

}

/**********************************************************************
 * Routine: ViscousFluxDiamondPath_n                                  *
 *                                                                    *
 * This routine computes the viscous flux at the specified quadrature *
 * point, X.  The gradient of the primitive variables is computed on  *
 * the diamond-path defined by the points X1, X2, X3, and X4.  Only   *
 * half of the diamond (three points) is used at solid or transpiring *
 * boundaries.                                                        *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState ViscousFluxDiamondPath_n(const Vector2D &X,
					const Vector2D &Xl, const Dusty2D_pState &Wl,
					const Vector2D &Xd, const Dusty2D_pState &Wd,
					const Vector2D &Xr, const Dusty2D_pState &Wr,
					const Vector2D &Xu, const Dusty2D_pState &Wu,
					const Vector2D &norm_dir,
					const int &Axisymmetric,
					const int &stencil_flag) {

  if (stencil_flag == DIAMONDPATH_NONE) return Dusty2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);

  Dusty2D_pState W_face, dWdxl, dWdyl, dWdxr, dWdyr, dWdx, dWdy;
  Dusty2D_cState Flux;
  double A, Al, Ar;
  Vector2D ndl, nud, nlu, nrd, nur, ndu;
  int error_flag;

  // Compute Green-Gauss integration on 'left' triangle:
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX ||
      stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
    ndl = Vector2D((Xd.y-Xl.y),-(Xd.x-Xl.x));
    nud = Vector2D((Xu.y-Xd.y),-(Xu.x-Xd.x));
    nlu = Vector2D((Xl.y-Xu.y),-(Xl.x-Xu.x));
    Al = HALF*((Xd-Xl)^(Xu-Xl));
    W_face = HALF*(Wd+Wl);
    dWdxl = W_face*ndl.x;
    dWdyl = W_face*ndl.y;
    W_face = HALF*(Wu+Wd);
    dWdxl += W_face*nud.x;
    dWdyl += W_face*nud.y;
    W_face = HALF*(Wl+Wu);
    dWdxl += W_face*nlu.x;
    dWdyl += W_face*nlu.y;
    dWdxl /= Al;
    dWdyl /= Al;
  }

  // Compute Green-Gauss integration on 'right' triangle:
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION ||
      stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX ||
      stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
    nrd = Vector2D((Xr.y-Xd.y),-(Xr.x-Xd.x));
    nur = Vector2D((Xu.y-Xr.y),-(Xu.x-Xr.x));
    ndu = Vector2D((Xd.y-Xu.y),-(Xd.x-Xu.x));
    Ar = HALF*((Xr-Xu)^(Xr-Xd));
    W_face = HALF*(Wr+Wd);
    dWdxr = W_face*nrd.x;
    dWdyr = W_face*nrd.y;
    W_face = HALF*(Wu+Wr);
    dWdxr += W_face*nur.x;
    dWdyr += W_face*nur.y;
    W_face = HALF*(Wd+Wu);
    dWdxr += W_face*ndu.x;
    dWdyr += W_face*ndu.y;
    dWdxr /= Ar;
    dWdyr /= Ar;
  }

  // Compute the viscous flux.
  if (stencil_flag == DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION) {
    // Determine the total area.
    A = Al + Ar;
    // Determine the area-averaged gradients of the primitive variables.
    dWdx = (Al*dWdxl + Ar*dWdxr)/A;
    dWdy = (Al*dWdyl + Ar*dWdyr)/A;
    // Determine the primitive state at the quadrature point.
    //W_face = (Wl + Wd + Wr + Wu)/FOUR;
    error_flag = Bilinear_Interpolation_ZY(Wl,Xl,Wu,Xu,Wr,Xr,Wd,Xd,X,W_face);
    //if (error_flag) return error_flag;
    Flux = ViscousFlux_n(X,W_face,dWdx,dWdy,norm_dir,Axisymmetric,OFF);

  } else if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX ||
	     stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
    W_face = Wd;
    if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX) {
      Flux = ViscousFlux_n(X,W_face,dWdxl,dWdyl,norm_dir,Axisymmetric,ON);
    } else if (stencil_flag == DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL) {
      Flux = ViscousFlux_n(X,W_face,dWdxl,dWdyl,norm_dir,Axisymmetric,OFF);
    }

  } else if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX ||
	     stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
    W_face = Wd;
    if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX) {
      Flux = ViscousFlux_n(X,W_face,dWdxr,dWdyr,norm_dir,Axisymmetric,ON);
    } else if (stencil_flag == DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL) {
      Flux = ViscousFlux_n(X,W_face,dWdxr,dWdyr,norm_dir,Axisymmetric,OFF);
    }

  }

  // Return the viscous flux.
  return Flux;

}

/**********************************************************************
 * Routine: ViscousFluxHybrid_n                                       *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 *                                                                    *
 **********************************************************************/
Dusty2D_cState ViscousFluxHybrid_n(const Vector2D &X,
				   Dusty2D_pState &W,
				   const Vector2D &X1,
				   const Dusty2D_pState &W1,
				   const Dusty2D_pState &dW1dx,
				   const Dusty2D_pState &dW1dy,
				   const Vector2D &X2,
				   const Dusty2D_pState &W2,
				   const Dusty2D_pState &dW2dx,
				   const Dusty2D_pState &dW2dy,
				   const Vector2D &norm_dir,
				   const int &Axisymmetric,
				   const int &adiabatic_flag) {

  Dusty2D_pState dWdx_ave, dWdy_ave, dWdx, dWdy, dWds;
  Vector2D dX;
  double ds;

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  dWdx_ave = HALF*(dW1dx + dW2dx);
  dWdy_ave = HALF*(dW1dy + dW2dy);

  dX = X2-X1; ds = dX.abs(); dX /= ds;

  dWds = (W2-W1)/ds;

  dWdx = dWdx_ave + (dWds - dWdx_ave*dX.x)*norm_dir.x/dot(norm_dir,dX);
  dWdy = dWdy_ave + (dWds - dWdy_ave*dX.y)*norm_dir.y/dot(norm_dir,dX);

  // Return the intermediate state solution viscous flux.
  return ViscousFlux_n(X,W,dWdx,dWdy,norm_dir,Axisymmetric,OFF);
  //return Dusty2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);

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
double WallShearStress(const Dusty2D_pState &W1,
		       const Vector2D &X1,
		       const Vector2D &X2,
		       const Vector2D &X3,
		       const Vector2D &nhat) {

  double A;
  Vector2D n21, n32, n13, that;
  Dusty2D_pState W2, W3, W_face, dWdx, dWdy, dWdt, dWdn;

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
  dWdn.v.x = dWdx.v.x*nhat.x + dWdy.v.x*nhat.y;
  dWdt.v.y = dWdx.v.y*that.x + dWdy.v.y*that.y;
  dWdn.v.y = dWdx.v.y*nhat.x + dWdy.v.y*nhat.y;

  // Return the wall shear stress.
  return W1.mu()*(dWdt.v.y + dWdn.v.x);

}

/**********************************************************************
 * Routine: WallShearStress                                           *
 *                                                                    *
 * This routine returns shear stress in the coordinate frame defined  *
 * by the given normal.                                               *
 *                                                                    *
 **********************************************************************/
double ShearStress(const Dusty2D_pState &W,
		   const Dusty2D_pState &dWdx,
		   const Dusty2D_pState &dWdy,
		   const Vector2D &nhat) {

  Vector2D that;
  Dusty2D_pState dWdt, dWdn;

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

/**********************************************************************
 * Routine: MultiVelocity_Switch                                      *
 **********************************************************************/
void Dusty2D_pState::MultiVelocity_Switch(void) {
  assert(NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION);
  Particle2D_pComponents Wcomp;
  // Consider each component of the particle family and switch
  // compoments if necessary.
  for (int nc = 0; nc < NUM_CMP_PART; nc++) {
    if (Wp[nc].sigma >= NANO) {
      if (fabs(Wp[nc].u.x) < NANO) Wp[nc].u.x = ZERO;
      if (fabs(Wp[nc].u.y) < NANO) Wp[nc].u.y = ZERO;
      if (Wp[nc].u.x > ZERO && Wp[nc].u.y >= ZERO) {
	Wcomp[0] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
      } else if (Wp[nc].u.x <= ZERO && Wp[nc].u.y > ZERO) {
	Wcomp[1] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
      } else if (Wp[nc].u.x < ZERO && Wp[nc].u.y <= ZERO) {
	Wcomp[2] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
      } else if (Wp[nc].u.x >= ZERO && Wp[nc].u.y < ZERO) {
	Wcomp[3] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
      }
    }
  }
  // Divide each (mass) state variable by the accumulated density.
  for (int nc = 0; nc < NUM_CMP_PART; nc++) {
    if (Wcomp[nc].sigma < NANO) {
      Wcomp[nc].Vacuum();
    } else {
      Wcomp[nc].u.x = Wcomp[nc].u.x/Wcomp[nc].sigma;
      Wcomp[nc].u.y = Wcomp[nc].u.y/Wcomp[nc].sigma;
      Wcomp[nc].Tp = Wcomp[nc].Tp/Wcomp[nc].sigma;
    }
  }
  // Copy the reset component state.
  Wp = Wcomp;
}

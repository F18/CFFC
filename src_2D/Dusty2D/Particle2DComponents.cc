/**********************************************************************
 * Particle2DComponents.cc: Subroutines for 2D particle components    *
 *                          solution classes.                         *
 **********************************************************************/

// Include 2D particle components solution header file.

#ifndef _PARTICLE2D_COMPONENTS_INCLUDED
#include "Particle2DComponents.h"
#endif // _PARTICLE2D_COMPONENTS_INCLUDED

/**********************************************************************
 * Particle2D_pComponents -- Create storage and assign default value. *
 **********************************************************************/
int Particle2D_pComponents::NUM_CMP_PARTICLE2D = PARTICLE2D_SINGLE_VELOCITY_FORMULATION;

/**********************************************************************
 * Particle2D_cComponents -- Create storage and assign default value. *
 **********************************************************************/
int Particle2D_cComponents::NUM_CMP_PARTICLE2D = PARTICLE2D_SINGLE_VELOCITY_FORMULATION;

/**********************************************************************
 * Routine: Rotate                                                    *
 *                                                                    *
 * This function returns the solution in the local rotated frame.     *
 *                                                                    *
 **********************************************************************/
Particle2D_pComponents Rotate(const Particle2D_pComponents &W,
			      const Vector2D &norm_dir) {

  Particle2D_pComponents Wrotated;
  Particle2D_pState Wtemp;

  // Determine the rotated components.
  if (W.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Wrotated[0].Rotate(W[0],norm_dir);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
      if (W[nc].sigma < NANO) {
      } else {
	Wtemp.Rotate(W[nc],norm_dir);
	Wrotated += Particle2D_pComponents(Wtemp);
      }
    }
  }

  // Return the rotated components.
  return Wrotated;

}

Particle2D_cComponents Rotate(const Particle2D_cComponents &U,
			      const Vector2D &norm_dir) {

  Particle2D_cComponents Urotated;
  Particle2D_cState Utemp;

  // Determine the rotated components.
  if (U.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Urotated[0].Rotate(U[0],norm_dir);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < U.NUM_CMP_PARTICLE2D; nc++) {
      Utemp.Rotate(U[nc],norm_dir);
      Urotated += Particle2D_cComponents(Utemp);
    }
  }

  // Return the rotated components.
  return Urotated;

}

/**********************************************************************
 * Routine: Translate                                                 *
 *                                                                    *
 * This function returns the solution in a stationary frame.          *
 *                                                                    *
 **********************************************************************/
Particle2D_pComponents Translate(const Particle2D_pComponents &W,
				 const Vector2D &V) {

  Particle2D_pComponents Wtranslated;
  Particle2D_pState Wtemp;

  // Determine the translated components.
  if (W.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Wtranslated[0].Translate(W[0],V);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
      if (W[nc].sigma < NANO) {
      } else {
	Wtemp.Translate(W[nc],V);
	Wtranslated += Particle2D_pComponents(Wtemp);
      }
    }
  }

  // Return the translated components.
  return Wtranslated;

}

/**********************************************************************
 * Routine: Reflect                                                   *
 *                                                                    *
 * This function returns the reflected solution state in a given      *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
Particle2D_pComponents Reflect(const Particle2D_pComponents &W,
			       const Vector2D &norm_dir) {

  Particle2D_pComponents Wreflected;
  Particle2D_pState Wtemp;

  // Determine the reflected components.
  if (W.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Wreflected[0].Reflect(W[0],norm_dir);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
      if (W[nc].sigma < NANO) {
      } else {
	//if (dot(norm_dir,W[nc].u) > ZERO) Wtemp.Reflect(W[nc],norm_dir);
	//else Wtemp.Copy(W[nc]);
	//Wreflected += Particle2D_pComponents(Wtemp);
	if (dot(norm_dir,W[nc].u) > ZERO) {
	  Wtemp.Reflect(W[nc],norm_dir);
	  Wreflected += Particle2D_pComponents(Wtemp);
	}
      }
    }
  }

  // Return the reflected components.
  return Wreflected;

}

/**********************************************************************
 * Routine: Mirror                                                    *
 *                                                                    *
 * This function returns the mirrored solution state in a given       *
 * direction given the primitive solution variables and the unit      *
 * normal vector in the direction of interest.                        *
 *                                                                    *
 **********************************************************************/
Particle2D_pComponents Mirror(const Particle2D_pComponents &W,
			      const Vector2D &norm_dir) {
  
  Particle2D_pComponents Wmirrored;
  Particle2D_pState Wtemp;

  // Determine the mirrored components.
  if (W.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Wmirrored[0].Mirror(W[0],norm_dir);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
      if (W[nc].sigma < NANO) {
      } else {
	Wtemp.Mirror(W[nc],norm_dir);
	Wmirrored += Particle2D_pComponents(Wtemp);
      }
    }
  }

  // Return the mirrored components.
  return Wmirrored;
  
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
Particle2D_pComponents Absorb(const Particle2D_pComponents &W,
			      const Vector2D &norm_dir) {
  
  Particle2D_pComponents Wabsorbed;
  Particle2D_pState Wtemp;

  // Determine the absorbed components.
  if (W.NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    Wabsorbed[0].Absorb(W[0],norm_dir);
  } else {
    // Multi-velocity formulation.
    for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
      if (W[nc].sigma < NANO) {
      } else {
	Wtemp.Absorb(W[nc],norm_dir);
	Wabsorbed += Particle2D_pComponents(Wtemp);
      }
    }
  }

  // Return the absorbed components.
  return Wabsorbed;

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
Particle2D_cComponents FluxSaurel_n(const Particle2D_pComponents &Wl,
				    const Particle2D_pComponents &Wr,
				    const Vector2D &norm_dir,
				    const double &cm) {

  Particle2D_cComponents Flux;

  // Determine the intermediate state solution flux.  It is assumed
  // that only one component exists.
  Flux[0].FluxSaurel_n(Wl[0],Wr[0],norm_dir,cm);

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
Particle2D_cComponents FluxSaurel_MB_n(const Particle2D_pComponents &Wl,
				       const Particle2D_pComponents &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir,
				       const double &cm) {

  Particle2D_cComponents Flux;

  // Determine the intermediate state solution flux.  It is assumed
  // that only one component exists.
  Flux[0].FluxSaurel_MB_n(Wl[0],Wr[0],V,norm_dir,cm);

  // Return the solution flux.
  return Flux;

}

/**********************************************************************
 * Routine: FluxMultiVelocity                                         *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by computing the  *
 * flux from the average of the left and right solution states.       *
 *                                                                    *
 **********************************************************************/
Particle2D_pComponents FluxMultiVelocity(const Particle2D_pComponents &Wl,
					 const Particle2D_pComponents &Wr) {

  Particle2D_pComponents W;

  // Determine the particle-phase flux for each component of each 
  // particle family.
  for (int nc = 0; nc < W.NUM_CMP_PARTICLE2D; nc++) {
    // Determine the intermediate state solution flux.
    if (Wl[nc].sigma < NANO && Wr[nc].sigma < NANO) {
      W[nc].Vacuum();
    } else {
      if (Wl[nc].u.x <= Wr[nc].u.x) {
	if ((Wl[nc].u.x >= ZERO) && (Wr[nc].u.x >= ZERO)) {
	  W[nc] = Wl[nc];
	} else if ((Wl[nc].u.x < ZERO) && (Wr[nc].u.x <= ZERO)) {
	  W[nc] = Wr[nc];
	} else {
	  cout << "ERROR MB 1" << endl; cout.flush();
	}
      } else if (Wl[nc].u.x > Wr[nc].u.x) {
	if ((Wl[nc].u.x > ZERO) && (Wr[nc].u.x >= ZERO)) {
	  W[nc] = Wl[nc];
	} else if ((Wl[nc].u.x <= ZERO) && (Wr[nc].u.x < ZERO)) {
	  W[nc] = Wr[nc];
	} else {
	  cout << "ERROR MB 2" << endl; cout.flush();
	}
      }
    }
    // Error checking.
    if (W[nc].sigma < NANO || W[nc].Tp < NANO) {
      W[nc].Vacuum();
    }
  }

  // Return the intermediate solution state.
  return W;
  
}

Particle2D_pComponents FluxMultiVelocity(const Particle2D_cComponents &Ul,
					 const Particle2D_cComponents &Ur,
					 const double &cm) {
  return FluxMultiVelocity(Ul.W(cm),Ur.W(cm));
}

/**********************************************************************
 * Routine: FluxMultiVelocity_n                                       *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by computing the  *
 * flux from the average of the left and right solution states.       *
 *                                                                    *
 **********************************************************************/
Particle2D_cComponents FluxMultiVelocity_n(const Particle2D_pComponents &Wl,
					   const Particle2D_pComponents &Wr,
					   const Vector2D &norm_dir,
					   const double &cm) {

//   Particle2D_pComponents Wl_rotated, Wr_rotated, W_rotated, W;
//   Particle2D_cComponents Flux, Flux_rotated;

//   // Apply the frame rotation and evaluate left and right solution 
//   // states in the local rotated frame defined by the unit normal 
//   // vector.
//   Wl_rotated.Rotate(Wl,norm_dir);
//   Wr_rotated.Rotate(Wr,norm_dir);

//   // Evaluate the intermediate solution state in the rotated frame.
//   W_rotated = FluxMultiVelocity(Wl_rotated,Wr_rotated);

//   // Rotate back to the original Cartesian reference frame
//   W.Rotate(W_rotated,Vector2D(norm_dir.x,-norm_dir.y));

//   // Determine the intermediate state solution flux.
//   Flux = W.Fn(cm,norm_dir);

//   Particle2D_pState Wl_rotated, Wr_rotated, W_rotated;
//   Particle2D_cState Flux_rotated;
//   Particle2D_cComponents Flux;

//   Wl_rotated.Rotate(Wl[0],norm_dir);
//   Wr_rotated.Rotate(Wr[0],norm_dir);
//   W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
//   Flux_rotated.F(W_rotated,cm);
//   Flux[0].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

//   Wl_rotated.Rotate(Wl[1],norm_dir);
//   Wr_rotated.Rotate(Wr[1],norm_dir);
//   W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
//   Flux_rotated.F(W_rotated,cm);
//   Flux[1].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

//   Wl_rotated.Rotate(Wl[2],norm_dir);
//   Wr_rotated.Rotate(Wr[2],norm_dir);
//   W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
//   Flux_rotated.F(W_rotated,cm);
//   Flux[2].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

//   Wl_rotated.Rotate(Wl[3],norm_dir);
//   Wr_rotated.Rotate(Wr[3],norm_dir);
//   W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
//   Flux_rotated.F(W_rotated,cm);
//   Flux[3].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  Particle2D_pState Wl_rotated, Wr_rotated, W_rotated;
  Particle2D_cState Flux_rotated;
  Particle2D_cComponents Flux;

  Wl_rotated.Rotate(Wl[0],norm_dir);
  Wr_rotated.Rotate(Wr[0],norm_dir);
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
  Flux_rotated.F(W_rotated,cm);
  Flux[0].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  Wl_rotated.Rotate(Wl[1],norm_dir);
  Wr_rotated.Rotate(Wr[1],norm_dir);
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
  Flux_rotated.F(W_rotated,cm);
  Flux[1].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  Wl_rotated.Rotate(Wl[2],norm_dir);
  Wr_rotated.Rotate(Wr[2],norm_dir);
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
  Flux_rotated.F(W_rotated,cm);
  Flux[2].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  Wl_rotated.Rotate(Wl[3],norm_dir);
  Wr_rotated.Rotate(Wr[3],norm_dir);
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);
  Flux_rotated.F(W_rotated,cm);
  Flux[3].Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the intermediate state solution flux.
  return Flux;

}

Particle2D_cComponents FluxMultiVelocity_n(const Particle2D_cComponents &Ul,
					   const Particle2D_cComponents &Ur,
					   const Vector2D &norm_dir,
					   const double &cm) {
  return FluxMultiVelocity_n(Ul.W(cm),Ur.W(cm),norm_dir,cm);
}

/**********************************************************************
 * Routine: FluxMultiVelocity_MB_n                                    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by computing the  *
 * flux from the average of the left and right solution states.       *
 *                                                                    *
 **********************************************************************/
Particle2D_cComponents FluxMultiVelocity_MB_n(const Particle2D_pComponents &Wl,
					      const Particle2D_pComponents &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir,
					      const double &cm) {

  Particle2D_pComponents Wl_rotated, Wr_rotated, W_rotated;
  Particle2D_cComponents Flux, Flux_rotated;
  Vector2D V_rotated;

  cout << endl << " ERROR: FluxMultiVelocity_MB_n not implemented.";
  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
//   Wl_rotated.Rotate(Wl,norm_dir);
//   Wr_rotated.Rotate(Wr,norm_dir);
//   V_rotated.Rotate(V,norm_dir);

//   // Evaluate the intermediate solution state in the rotated frame.
//   W_rotated = FluxMultiVelocity(Wl_rotated,Wr_rotated);

//   // Evaluate the intermediate state solution flux in the rotated 
//   // frame.
//   Flux_rotated = W_rotated.F(cm);

//   // Rotate back to the original Cartesian reference frame and return 
//   // the solution flux.
//   Flux.Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

  // Return the solution flux.
  return Flux;

}

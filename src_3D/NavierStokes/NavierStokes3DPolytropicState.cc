/*! \file  NavierStokes3DPolytropicState.cc
 *	\brief Subroutines for 3D Navier Stokes Solution State Classes. 
 */

/* Include 3D NavierStokes solution state header file. */

#include "NavierStokes3DPolytropicState.h"


/* -------------------------------------------------------------------------- *
 *				 NavierStokes3D_Polytropic_pState subroutines				  *
 * -------------------------------------------------------------------------- */


/*
 * Gas specific constants
 * ----------------------
 */
double NavierStokes3D_Polytropic_pState::Cp = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double NavierStokes3D_Polytropic_pState::Cv = (R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double NavierStokes3D_Polytropic_pState::v1 = AIR_c1; 
double NavierStokes3D_Polytropic_pState::v2 = AIR_c2; 
double NavierStokes3D_Polytropic_pState::v3 = AIR_c3;
double NavierStokes3D_Polytropic_pState::v4 = AIR_c4; 
double NavierStokes3D_Polytropic_pState::v5 = AIR_c5;

/* 
 * Set gas specific constants
 * --------------------------
 */
/*!
 * Sets gas type and corresponding static variables to "AIR"
 */
void NavierStokes3D_Polytropic_pState::setgas(void) {
	gas_type = "AIR";
	g  = GAMMA_AIR;
	R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
	v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
	gm1 = g-ONE;
	gm1i = ONE/gm1;
	Cp = g*R*gm1i;
	Cv = R*gm1i;
}
/*!
 * Sets gas type and corresponding static variables to the specified gas
 * \param[in] string_ptr name of the gas (e.g. "AIR", "H2", "HE", "N2", "O2")
 */		
void NavierStokes3D_Polytropic_pState::setgas(char *str_ptr) {
    gas_type = str_ptr;
	if (strcmp(str_ptr,"AIR") == 0) {
		g  = GAMMA_AIR;
		R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
	} else if (strcmp(str_ptr,"H2") == 0) {
		g  = GAMMA_H2;
		R  = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
		v1 = H2_c1; v2 = H2_c2; v3 = H2_c3; v4 = H2_c4; v5 = H2_c5;
	} else if (strcmp(str_ptr,"HE") == 0) {
		g  = GAMMA_HE;
		R  = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
		v1 = HE_c1; v2 = HE_c2; v3 = HE_c3; v4 = HE_c4; v5 = HE_c5;
	} else if (strcmp(str_ptr,"N2") == 0) {
		g  = GAMMA_N2;
		R  = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
		v1 = N2_c1; v2 = N2_c2; v3 = N2_c3; v4 = N2_c4; v5 = N2_c5;
	} else if (strcmp(str_ptr,"O2") == 0) {
		g  = GAMMA_O2;
		R  = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
		v1 = O2_c1; v2 = O2_c2; v3 = O2_c3; v4 = O2_c4; v5 = O2_c5;
	} else {
		gas_type = "AIR";
		g  = GAMMA_AIR;
		R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; AIR_c4;
		cerr << "\n Gas type '"<< gas_type << "' not supported for"
			 << "Navier Stokes 3D in NavierStokes3DPolytropicState.cc \n"
			 << "Gas type set to 'AIR'." << endl;
	}
	gm1 = g-ONE;
	gm1i = ONE/gm1;
	Cp = g*R*gm1i;
	Cv = R*gm1i;
}


/*
 * Viscous- and Heat Fluxes
 * ------------------------
 */
// x- direction viscous flux and heat flux
/*! \f[
        Gx = \left( \begin{array}{c}
                    0 \\ \tau_{xx} \\ \tau_{xy} \\ \tau_{xz} \\
                    u \tau_{xx} + v \tau_{xy} + w \tau_{xz} \end{array} \right)
           + \left( \begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ -q_x \end{array} \right)
    \f]
*/
NavierStokes3D_Polytropic_cState NavierStokes3D_Polytropic_pState::Fvx(const NavierStokes3D_Polytropic_pState &dWdx,
                                                                       const NavierStokes3D_Polytropic_pState &dWdy,
                                                                       const NavierStokes3D_Polytropic_pState &dWdz) {	
	Tensor3D tau = tau_x(dWdx,dWdy,dWdz);
	Vector3D q = q_x(dWdx,dWdy,dWdz);
	return (NavierStokes3D_Polytropic_cState(ZERO,tau.xx,tau.xy,tau.xz,v.x*tau.xx+v.y*tau.xy+v.z*tau.xz-q.x));
}

NavierStokes3D_Polytropic_cState NavierStokes3D_Polytropic_pState::Fv(const NavierStokes3D_Polytropic_pState &dWdx,
                                                                      const NavierStokes3D_Polytropic_pState &dWdy,
                                                                      const NavierStokes3D_Polytropic_pState &dWdz) {	
	Tensor3D tau = tau_x(dWdx,dWdy,dWdz);
	Vector3D q = q_x(dWdx,dWdy,dWdz);
	return (NavierStokes3D_Polytropic_cState(ZERO,tau.xx,tau.xy,tau.xz,v.x*tau.xx+v.y*tau.xy+v.z*tau.xz-q.x));
}


// y-direction viscous flux and heat flux
/*! \f[
Gy = \left( \begin{array}{c}
            0 \\ \tau_{xy} \\ \tau_{yy} \\ \tau_{yz} \\
            u \tau_{xy} + v \tau_{yy} + w \tau_{yz} \end{array} \right)
+ \left( \begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ -q_y \end{array} \right)
    \f]
    */
NavierStokes3D_Polytropic_cState NavierStokes3D_Polytropic_pState::Fvy(const NavierStokes3D_Polytropic_pState &dWdx,
																	  const NavierStokes3D_Polytropic_pState &dWdy,
																	  const NavierStokes3D_Polytropic_pState &dWdz) {	
	Tensor3D tau = tau_y(dWdx,dWdy,dWdz);
	Vector3D q = q_y(dWdx,dWdy,dWdz);
	return (NavierStokes3D_Polytropic_cState(ZERO,tau.xy,tau.yy,tau.yz,v.x*tau.xy+v.y*tau.yy+v.z*tau.yz-q.y));
}


// z-direction viscous flux and heat flux
/*! \f[
Gz = \left( \begin{array}{c}
            0 \\ \tau_{xz} \\ \tau_{yz} \\ \tau_{zz} \\
            u \tau_{xz} + v \tau_{yz} + w \tau_{zz} \end{array} \right)
+ \left( \begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ -q_z \end{array} \right)
    \f]
    */
NavierStokes3D_Polytropic_cState NavierStokes3D_Polytropic_pState::Fvz(const NavierStokes3D_Polytropic_pState &dWdx,
																	  const NavierStokes3D_Polytropic_pState &dWdy,
																	  const NavierStokes3D_Polytropic_pState &dWdz) {	
	Tensor3D tau = tau_z(dWdx,dWdy,dWdz);
	Vector3D q = q_z(dWdx,dWdy,dWdz);
	return (NavierStokes3D_Polytropic_cState(ZERO,tau.xz,tau.yz,tau.zz,v.x*tau.xz+v.y*tau.yz+v.z*tau.zz-q.z));
}



// n-direction viscous flux and heat flux
NavierStokes3D_Polytropic_cState  NavierStokes3D_Polytropic_pState::FluxViscous_n(
			  const NavierStokes3D_Polytropic_pState &Wl,
			  const NavierStokes3D_Polytropic_pState &Wr,
			  const NavierStokes3D_Polytropic_pState &Wc,
			  const NavierStokes3D_Polytropic_pState &Wc_Neighbour,
			  const NavierStokes3D_Polytropic_pState &dWdx,
			  const NavierStokes3D_Polytropic_pState &dWdy,
			  const NavierStokes3D_Polytropic_pState &dWdz,
			  const NavierStokes3D_Polytropic_pState &dWdx_Neighbour,
			  const NavierStokes3D_Polytropic_pState &dWdy_Neighbour,
			  const NavierStokes3D_Polytropic_pState &dWdz_Neighbour,
			  const Vector3D &norm, const Vector3D &ts, const double &deltad,
			  const double &Volume, const double &Volume_Neighbour){

																				  	
	// construct the gradients on the cell interface (surface) 
	// based on Hybrid Average Gradient-Diamond-Path Approach
	// Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
	//            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]
	
	// weighted factor based on volume
	double alpha = Volume/(Volume + Volume_Neighbour);
	
	
	NavierStokes3D_Polytropic_pState dWdx_Weighted, 
		dWdy_Weighted, dWdz_Weighted, dWdx_face, 
		dWdy_face, dWdz_face, Grad_middle_term;
	
	NavierStokes3D_Polytropic_pState W_face;
	
	dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neighbour;
	dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neighbour;
	dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neighbour;
	
	
	// a weighted term  
	Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y +
		dWdz_Weighted*ts.z;
	
	// gradients of primitive variables on the face
	dWdx_face = (Wc_Neighbour - Wc)/deltad *norm.x/dot(norm, ts) + (dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
	dWdy_face = (Wc_Neighbour - Wc)/deltad *norm.y/dot(norm, ts) + (dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
	dWdz_face = (Wc_Neighbour - Wc)/deltad *norm.z/dot(norm, ts) + (dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));
	
	W_face = HALF*(Wl + Wr);
	
    if (fabs(norm.y) < TOLER && fabs(norm.z) < TOLER) {
        return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face)*norm.x);
        
    } else if (fabs(norm.x) < TOLER && fabs(norm.z) < TOLER) {
        return (W_face.Fvy(dWdx_face, dWdy_face, dWdz_face)*norm.y);
        
    } else if (fabs(norm.x) < TOLER && fabs(norm.y) < TOLER) {
        return (W_face.Fvz(dWdx_face, dWdy_face, dWdz_face)*norm.z);
        
    } else {
        return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face)*norm.x +
                W_face.Fvy(dWdx_face, dWdy_face, dWdz_face)*norm.y +
                W_face.Fvz(dWdx_face, dWdy_face, dWdz_face)*norm.z);
    }
}

/*
 * Navier-Stokes related functions
 * -------------------------------
 */

// Dynamic viscosity
/*!
 * Viscosity calculated using a curve fit depending on temperature
 * mu = mu_gottlieb(v1,v2,v3,v4,v5,T) defined in GasConstants.h
 */
double NavierStokes3D_Polytropic_pState::mu(void) {		
	return (mu_gottlieb(v1,v2,v3,v4,v5,T()));
}


// Kinematic viscosity
/*!
 * nu = mu()/rho 
 */
double NavierStokes3D_Polytropic_pState::nu(void) {		
	return (mu()/rho);
}

// thermal conductivity
/*!
 * Thermal conductivity calculated using a curve fit depending on temperature
 * k = kappa_gottlieb(v1,v2,v3,v4,v5,T,g,Cp) defined in GasConstants.h
 */
double NavierStokes3D_Polytropic_pState::kappa(void) {
	return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,Cp);
}

// Molecular stress
Tensor3D NavierStokes3D_Polytropic_pState::tau(const NavierStokes3D_Polytropic_pState &dWdx, 
                                               const NavierStokes3D_Polytropic_pState &dWdy,
                                               const NavierStokes3D_Polytropic_pState &dWdz) {
	Tensor3D tau;
    double mu_=mu();
	tau.xx = 1.0/3.0 *mu_*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z);
	tau.xy = mu_*(dWdx.v.y + dWdy.v.x);
	tau.xz = mu_*(dWdx.v.z + dWdz.v.x);
	tau.yy = 1.0/3.0 *mu()*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z);
	tau.yz = mu_*(dWdy.v.z + dWdz.v.y);
	tau.zz = 1.0/3.0 *mu_*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y);
    
	return (tau);	
}

Tensor3D NavierStokes3D_Polytropic_pState::tau_x(const NavierStokes3D_Polytropic_pState &dWdx, 
                                                 const NavierStokes3D_Polytropic_pState &dWdy,
                                                 const NavierStokes3D_Polytropic_pState &dWdz) {
	Tensor3D tau;
    double mu_=mu();
	tau.xx = 1.0/3.0 *mu_*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z);
	tau.xy = mu_*(dWdx.v.y + dWdy.v.x);
	tau.xz = mu_*(dWdx.v.z + dWdz.v.x);    
	return (tau);	
}

Tensor3D NavierStokes3D_Polytropic_pState::tau_y(const NavierStokes3D_Polytropic_pState &dWdx, 
                                                 const NavierStokes3D_Polytropic_pState &dWdy,
                                                 const NavierStokes3D_Polytropic_pState &dWdz) {
	Tensor3D tau;
    double mu_=mu();
	tau.xy = mu_*(dWdx.v.y + dWdy.v.x);
	tau.yy = 1.0/3.0 *mu()*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z);
	tau.yz = mu_*(dWdy.v.z + dWdz.v.y);
    
	return (tau);	
}

Tensor3D NavierStokes3D_Polytropic_pState::tau_z(const NavierStokes3D_Polytropic_pState &dWdx, 
                                                 const NavierStokes3D_Polytropic_pState &dWdy,
                                                 const NavierStokes3D_Polytropic_pState &dWdz) {
	Tensor3D tau;
    double mu_=mu();
	tau.xz = mu_*(dWdx.v.z + dWdz.v.x);
	tau.yz = mu_*(dWdy.v.z + dWdz.v.y);
	tau.zz = 1.0/3.0 *mu_*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y);
    
	return (tau);	
}





// Heat flux vector
/*!
 * q.x = -k dT/dx
 * q.y = -k dT/dy
 * q.z = -k dT/dz
 */
Vector3D NavierStokes3D_Polytropic_pState::q(const NavierStokes3D_Polytropic_pState &dWdx, 
                                             const NavierStokes3D_Polytropic_pState &dWdy,
                                             const NavierStokes3D_Polytropic_pState &dWdz) {
	double kappa_ = kappa();
	Vector3D q;
	q.x = -kappa_*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
	q.y = -kappa_*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);
	q.z = -kappa_*(dWdz.p - (p/rho)*dWdz.rho)/(rho*R);
	return (q);
}

Vector3D NavierStokes3D_Polytropic_pState::q_x(const NavierStokes3D_Polytropic_pState &dWdx, 
                                               const NavierStokes3D_Polytropic_pState &dWdy,
                                               const NavierStokes3D_Polytropic_pState &dWdz) {
	Vector3D q;
	q.x = -kappa()*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
	return (q);
}

Vector3D NavierStokes3D_Polytropic_pState::q_y(const NavierStokes3D_Polytropic_pState &dWdx, 
                                               const NavierStokes3D_Polytropic_pState &dWdy,
                                               const NavierStokes3D_Polytropic_pState &dWdz) {
	Vector3D q;
	q.y = -kappa()*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);
	return (q);
}

Vector3D NavierStokes3D_Polytropic_pState::q_z(const NavierStokes3D_Polytropic_pState &dWdx, 
                                               const NavierStokes3D_Polytropic_pState &dWdy,
                                               const NavierStokes3D_Polytropic_pState &dWdz) {
	Vector3D q;
	q.z = -kappa()*(dWdz.p - (p/rho)*dWdz.rho)/(rho*R);
	return (q);
}



/* -------------------------------------------------------------------------- *
 *				 NavierStokes3D_Polytropic_cState subroutines				  *
 * -------------------------------------------------------------------------- */


/*
 * Gas specific constants
 * ----------------------
 */
double NavierStokes3D_Polytropic_cState::Cp = (GAMMA_AIR*R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double NavierStokes3D_Polytropic_cState::Cv = (R_UNIVERSAL/(MOLE_WT_AIR*MILLI))/(GAMMA_AIR-ONE);
double NavierStokes3D_Polytropic_cState::v1 = AIR_c1; 
double NavierStokes3D_Polytropic_cState::v2 = AIR_c2; 
double NavierStokes3D_Polytropic_cState::v3 = AIR_c3;
double NavierStokes3D_Polytropic_cState::v4 = AIR_c4; 
double NavierStokes3D_Polytropic_cState::v5 = AIR_c5;


/* 
 * Set gas specific constants
 * --------------------------
 */
/*!
 * Sets gas type and corresponding static variables to "AIR"
 */
void NavierStokes3D_Polytropic_cState::setgas(void) {
	gas_type = "AIR";
	g  = GAMMA_AIR;
	R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
	v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
	gm1 = g-ONE;
	gm1i = ONE/gm1;
	Cp = g*R*gm1i;
	Cv = R*gm1i;
}
/*!
 * Sets gas type and corresponding static variables to the specified gas
 * \param[in] string_ptr name of the gas (e.g. "AIR", "H2", "HE", "N2", "O2")
 */	
void NavierStokes3D_Polytropic_cState::setgas(char *str_ptr) {
    gas_type = str_ptr;
	if (strcmp(str_ptr,"AIR") == 0) {
		g  = GAMMA_AIR;
		R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
	} else if (strcmp(str_ptr,"H2") == 0) {
		g  = GAMMA_H2;
		R  = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
		v1 = H2_c1; v2 = H2_c2; v3 = H2_c3; v4 = H2_c4; v5 = H2_c5;
	} else if (strcmp(str_ptr,"HE") == 0) {
		g  = GAMMA_HE;
		R  = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
		v1 = HE_c1; v2 = HE_c2; v3 = HE_c3; v4 = HE_c4; v5 = HE_c5;
	} else if (strcmp(str_ptr,"N2") == 0) {
		g  = GAMMA_N2;
		R  = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
		v1 = N2_c1; v2 = N2_c2; v3 = N2_c3; v4 = N2_c4; v5 = N2_c5;
	} else if (strcmp(str_ptr,"O2") == 0) {
		g  = GAMMA_O2;
		R  = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
		v1 = O2_c1; v2 = O2_c2; v3 = O2_c3; v4 = O2_c4; v5 = O2_c5;
	} else {
		gas_type = "AIR";
		g  = GAMMA_AIR;
		R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; AIR_c4;
		cerr << "\n Gas type '"<< gas_type << "' not supported for"
            << "Navier Stokes 3D in NavierStokes3DPolytropicState.cc \n"
            << "Gas type set to 'AIR'." << endl;
	}
	gm1 = g-ONE;
	gm1i = ONE/gm1;
	Cp = g*R*gm1i;
	Cv = R*gm1i;
}




/* Euler3DPolytropicState.cc:  Subroutines for 3D Euler Solution State Classes. */

/* Include 3D Euler solution state header file. */

#include "Euler3DPolytropicState.h"


/* -------------------------------------------------------------------------- *
 *				 Euler3D_Polytropic_pState subroutines						  *
 * -------------------------------------------------------------------------- */

   /*
	* Constructors
	* ------------
	*/

	/*! Creates primitive state with standard atmosphere variables */
	Euler3D_Polytropic_pState::Euler3D_Polytropic_pState() {
		rho = DENSITY_STDATM;	v.zero();	 p = PRESSURE_STDATM;
	}

	/*! Copies primitive state from a given primitive state */
	Euler3D_Polytropic_pState::Euler3D_Polytropic_pState(const Euler3D_Polytropic_pState &W) {
		rho = W.rho;	v = W.v;	p = W.p;
	}

	/*! Creates primitive state with given variables */ 
	Euler3D_Polytropic_pState::Euler3D_Polytropic_pState(const double &d, 
							  const Vector3D &V, 
							  const double &pre) {
		rho = d;	v = V;		p = pre;
	}

	/*! Creates primitive state with given variables */ 
	Euler3D_Polytropic_pState::Euler3D_Polytropic_pState(const double &d, 
							  const double &vx, const double &vy, const double &vz, 
							  const double &pre) {
		rho = d;	
		v.x = vx;	v.y = vy;	v.z = vz;
		p = pre;
	}

	/*! Creates primitive state from a given conservative state */ 
	Euler3D_Polytropic_pState::Euler3D_Polytropic_pState(const Euler3D_Polytropic_cState &U) {
		rho = U.rho;	v = U.v();	p = U.p();
	}

   /*
	* Set static variables 
	* -------------------- 
	*/
	
	int Euler3D_Polytropic_pState::num_vars = NUM_VAR_EULER3D;
	double Euler3D_Polytropic_pState::g = GAMMA_AIR;
	double Euler3D_Polytropic_pState::gm1 = GAMMA_AIR-ONE;
	double Euler3D_Polytropic_pState::gm1i = ONE/(GAMMA_AIR-ONE);
	double Euler3D_Polytropic_pState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);


	// Set gas constants.
	void Euler3D_Polytropic_pState::setgas(void) {
		g = GAMMA_AIR;
		R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		gm1 = g - ONE;
		gm1i = ONE/gm1;
	}
	void Euler3D_Polytropic_pState::setgas(char *string_ptr) {
		if (strcmp(string_ptr, "AIR") == 0) {
			g = GAMMA_AIR;
			R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		} else if (strcmp(string_ptr, "A") == 0) {
			g = GAMMA_A;
			R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
		} else if (strcmp(string_ptr, "CO") == 0) {
			g = GAMMA_CO;
			R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
		} else if (strcmp(string_ptr, "CO2") == 0) {
			g = GAMMA_CO2;
			R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
		} else if (strcmp(string_ptr, "CH4") == 0) {
			g = GAMMA_CH4;
			R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
		} else if (strcmp(string_ptr, "H") == 0) {
			g = GAMMA_H;
			R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
		} else if (strcmp(string_ptr, "H2") == 0) {
			g = GAMMA_H2;
			R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
		} else if (strcmp(string_ptr, "HE") == 0) {
			g = GAMMA_HE;
			R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
		} else if (strcmp(string_ptr, "H2O") == 0) {
			g = GAMMA_H2O;
			R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
		} else if (strcmp(string_ptr, "N2") == 0) {
			g = GAMMA_N2;
			R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
		} else if (strcmp(string_ptr, "O") == 0) {
			g = GAMMA_O;
			R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
		} else if (strcmp(string_ptr, "O2") == 0) {
			g = GAMMA_O2;
			R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
		} else if (strcmp(string_ptr, "e") == 0) {
			g = GAMMA_e;
			R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
		} else {
			g = GAMMA_AIR;
			R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		} 
	   gm1 = g - ONE;
	   gm1i = ONE/gm1;
	}

	// Total velocity.
	double Euler3D_Polytropic_pState::uo(void) {
		return (sqrt(v.x*v.x + v.y*v.y + v.z*v.z));  
	}
	double Euler3D_Polytropic_pState::uo(void) const{
		return (sqrt(v.x*v.x + v.y*v.y + v.z*v.z));  
	}

	// Temperature.
	double Euler3D_Polytropic_pState::T(void) {
		return (p/(rho*R));
	}
	double Euler3D_Polytropic_pState::T(void) const {
		return (p/(rho*R));
	}

	// gasconstant. (for compatibility reasons in hexablock.h)
	double Euler3D_Polytropic_pState::Rtot(void) {
		return (R);
	}
	double Euler3D_Polytropic_pState::Rtot(void) const {
		return (R);
	}
	// Specific internal energy.
	double Euler3D_Polytropic_pState::e(void) {
		return (p/(gm1*rho));
	}
	double Euler3D_Polytropic_pState::e(void) const {
		return (p/(gm1*rho));
	}

	// Total energy.
	double Euler3D_Polytropic_pState::E(void) {
		return (p*gm1i + HALF*rho*v.sqr());
	}
	double Euler3D_Polytropic_pState::E(void) const {
		return (p*gm1i + HALF*rho*v.sqr());
	}

	// Specific enthalpy. 
	double Euler3D_Polytropic_pState::h(void) {
		return (g*gm1i*p/rho + HALF*v.sqr());
	}
	double Euler3D_Polytropic_pState::h(void) const {
		return (g*gm1i*p/rho + HALF*v.sqr());
	}

	// Total enthalpy. 
	double Euler3D_Polytropic_pState::H(void) {
		return (g*gm1i*p + HALF*rho*v.sqr());
	}
	double Euler3D_Polytropic_pState::H(void) const {
		return (g*gm1i*p + HALF*rho*v.sqr());
	}

	// Sound speed. 
	double Euler3D_Polytropic_pState::a(void) {
		return (sqrt(g*p/rho));
	}
	double Euler3D_Polytropic_pState::a(void) const {
		return (sqrt(g*p/rho));
	}

	// Sound speed squared. 
	double Euler3D_Polytropic_pState::a2(void) {
		return (g*p/rho);
	}
	double Euler3D_Polytropic_pState::a2(void) const {
		return (g*p/rho);
	}

	// Mach number. 
	double Euler3D_Polytropic_pState::M(void) {
		return (uo()/sqrt(g*p/rho));
	}
	double Euler3D_Polytropic_pState::M(void) const {
		return (uo()/sqrt(g*p/rho));
	}

	// Specific entropy. 
	double Euler3D_Polytropic_pState::s(void) {
		return (R*gm1i*log(p/pow(rho, g)));
	}
	double Euler3D_Polytropic_pState::s(void) const {
		return (R*gm1i*log(p/pow(rho, g)));
	}

	// Momentum. 
	Vector3D Euler3D_Polytropic_pState::rhov(void) {
		return (rho*v);
	}
	Vector3D Euler3D_Polytropic_pState::rhov(void) const {
		return (rho*v);
		
	}
	double Euler3D_Polytropic_pState::rhov(const Vector3D &n) {
		return (rho*(v*n));
		
	}
	double Euler3D_Polytropic_pState::rhov(const Vector3D &n) const {
		return (rho*(v*n));
		
	}

	// Stagnation temperature. 
	double Euler3D_Polytropic_pState::To(void) {
		return ((p/(rho*R))*(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
	}
	double Euler3D_Polytropic_pState::To(void) const {
		return ((p/(rho*R))*(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
	}

	// Stagnation pressure. 
	double Euler3D_Polytropic_pState::po(void) {
		return (p*pow(ONE+HALF*gm1*sqr(uo())/(g*p/rho), g*gm1i));
	}
	double Euler3D_Polytropic_pState::po(void) const {
		return (p*pow(ONE+HALF*gm1*sqr(uo())/(g*p/rho), g*gm1i));
	}

	// Stagnation sound speed. 
	double Euler3D_Polytropic_pState::ao(void) {
		return (sqrt((g*p/rho)*(ONE+HALF*gm1*sqr(uo())/(g*p/rho))));
	}
	double Euler3D_Polytropic_pState::ao(void) const {
		return (sqrt((g*p/rho)*(ONE+HALF*gm1*sqr(uo())/(g*p/rho))));
	}

	// Stagnation enthalpy. 
	double Euler3D_Polytropic_pState::ho(void) {
		return ((g*p/(gm1*rho) + HALF*sqr(uo()))*
				(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
	}
	double Euler3D_Polytropic_pState::ho(void) const {
		return ((g*p/(gm1*rho) + HALF*sqr(uo()))
				*(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
	}

   /*
	* Operators. 
	* ---------- 
	*/
	
	// Index operator.
	double& Euler3D_Polytropic_pState::operator[](int index) {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (rho);
			case 2 :
				return (v.x);
			case 3 :
				return (v.y);
			case 4 :
				return (v.z);
			case 5 :
				return (p);
			default:
				return (rho);
		}
	}
	const double& Euler3D_Polytropic_pState::operator[](int index) const {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (rho);
			case 2 :
				return (v.x);
			case 3 :
				return (v.y);
			case 4 :
				return (v.z);
			case 5 :
				return (p);
			default:
				return (rho);
		}
	}

	// Binary arithmetic operators.
	Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho+W2.rho,W1.v+W2.v,W1.p+W2.p));
	}
	Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho-W2.rho,W1.v-W2.v, W1.p-W2.p));
	}
	double operator *(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho*W2.rho + W1.v.x*W2.v.x + W1.v.y*W2.v.y + W1.v.z*W2.v.z + W1.p*W2.p);
	}	
	Euler3D_Polytropic_pState operator *(const Euler3D_Polytropic_pState &W, const double &a) {
		return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z, a*W.p));
	}
	Euler3D_Polytropic_pState operator *(const double &a, const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z,a*W.p));
	}
	Euler3D_Polytropic_pState operator /(const Euler3D_Polytropic_pState &W, const double &a) {
		return (Euler3D_Polytropic_pState(W.rho/a, W.v.x/a, W.v.y/a, W.v.z/a,W.p/a));
	}
	Euler3D_Polytropic_pState operator ^(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho*W2.rho,
										  W1.v.x*W2.v.x, 
										  W1.v.y*W2.v.y, 
										  W1.v.z*W2.v.z, 
										  W1.p*W2.p));
	}
	
	// Unary arithmetic operators.
	Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_pState(W.rho,W.v.x,W.v.y, W.v.z,W.p));
	}
	Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_pState(-W.rho,-W.v.x, -W.v.y, -W.v.z, -W.p));
	}

	// Shortcut arithmetic operators.
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator +=(const Euler3D_Polytropic_pState &W) {
		rho += W.rho; 
		v.x += W.v.x;	v.y += W.v.y;	v.z += W.v.z;
		p += W.p;
		return *this;
	}
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator -=(const Euler3D_Polytropic_pState &W) {
		rho -= W.rho; 
		v.x -= W.v.x;	v.y -= W.v.y;	v.z -= W.v.z;
		p -= W.p;
		return *this;
	}
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator *=(const double &a) {
		rho *= a; 
		v.x *= a;	v.y *= a;	v.z *= a; 
		p *= a;
		return *this;
	}
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator /=(const double &a) {
		rho /= a; 
		v.x /= a;	v.y /= a;	v.z /= a;
		p /= a;
		return *this;
	}

	// Relational operators.
	int operator ==(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho == W2.rho && W1.v.x == W2.v.x && W1.v.y == W2.v.y
				&& W1.v.z == W2.v.z && W1.p == W2.p);
	}

	int operator !=(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho != W2.rho || W1.v.x != W2.v.x || W1.v.y != W2.v.y
				|| W1.v.z != W2.v.z || W1.p != W2.p);
	}

	// Input-output operators.
	ostream& operator << (ostream &out_file, const Euler3D_Polytropic_pState &W) {
		out_file.setf(ios::scientific);
		out_file << " " << W.rho<< " " << W.v.x << " " << W.v.y 
			<< " " << W.v.z << " " << W.p;
		out_file.unsetf(ios::scientific);
		return (out_file);
	}
	istream& operator >> (istream &in_file,  Euler3D_Polytropic_pState &W) {
		in_file.setf(ios::skipws);
		in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z 
			>> W.p;
		in_file.unsetf(ios::skipws);
		return (in_file);
	}


   /* 
	* Conserved solution state. 
	* -------------------------
	*/

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::U(void) {
		return (Euler3D_Polytropic_cState(rho, rhov(), E()));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::U(void) const {
		return (Euler3D_Polytropic_cState(rho, rhov(), E()));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::U(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho, W.rhov(), W.E()));
	}

   /* 
	* Fluxes and Jacobians
	* --------------------
	*/

	// x-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fx(void) {
		return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fx(void) const {
		return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
	}
	Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho*W.v.x, W.rho*sqr(W.v.x) + W.p, W.rho*W.v.x*W.v.y, W.rho*W.v.x*W.v.z, W.v.x*W.H()));
	}
	void Euler3D_Polytropic_pState::dFxdU(DenseMatrix &dFxdU) {
		dFxdU(0,1) += ONE;
		dFxdU(1,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.x));
		dFxdU(1,1) -= v.x*(g-THREE);
		dFxdU(1,2) -= v.y*gm1;
		dFxdU(1,3) -= v.z*gm1;
		dFxdU(1,4) += gm1;
		dFxdU(2,0) -= v.x*v.y;
		dFxdU(2,1) += v.y;
		dFxdU(2,2) += v.x;
		dFxdU(3,0) -= v.x*v.z;
		dFxdU(3,1) += v.z;
		dFxdU(3,3) += v.x;
		dFxdU(4,0) += HALF*v.x*v.sqr()*(g-TWO) - v.x*p/rho*g*gm1i;
		dFxdU(4,1) += HALF*(v.sqr() - TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
		dFxdU(4,2) -= v.x*v.y*gm1;
		dFxdU(4,3) -= v.x*v.z*gm1;
		dFxdU(4,4) += g*v.x;
	}
	void Euler3D_Polytropic_pState::dFxdU(DenseMatrix &dFxdU) const {
		dFxdU(0,1) += ONE;
		dFxdU(1,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.x));
		dFxdU(1,1) -= v.x*(g-THREE);
		dFxdU(1,2) -= v.y*gm1;
		dFxdU(1,3) -= v.z*gm1;
		dFxdU(1,4) += gm1;
		dFxdU(2,0) -= v.x*v.y;
		dFxdU(2,1) += v.y;
		dFxdU(2,2) += v.x;
		dFxdU(3,0) -= v.x*v.z;
		dFxdU(3,1) += v.z;
		dFxdU(3,3) += v.x;
		dFxdU(4,0) += HALF*v.x*v.sqr()*(g-TWO) - v.x*p/rho*g*gm1i;
		dFxdU(4,1) += HALF*(v.sqr() - TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
		dFxdU(4,2) -= v.x*v.y*gm1;
		dFxdU(4,3) -= v.x*v.z*gm1;
		dFxdU(4,4) += g*v.x;
	}
	void Euler3D_Polytropic_pState::dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_pState &W) {
		dFxdU(0,1) += ONE;
		dFxdU(1,0) += HALF*(W.v.sqr()*W.gm1 - TWO*sqr(W.v.x));
		dFxdU(1,1) -= W.v.x*(W.g-THREE);
		dFxdU(1,2) -= W.v.y*W.gm1;
		dFxdU(1,3) -= W.v.z*W.gm1;
		dFxdU(1,4) += W.gm1;
		dFxdU(2,0) -= W.v.x*W.v.y;
		dFxdU(2,1) += W.v.y;
		dFxdU(2,2) += W.v.x;
		dFxdU(3,0) -= W.v.x*W.v.z;
		dFxdU(3,1) += W.v.z;
		dFxdU(3,3) += W.v.x;
		dFxdU(4,0) += HALF*W.v.x*W.v.sqr()*(W.g-TWO) - W.v.x*W.p/W.rho*W.g*W.gm1i;
		dFxdU(4,1) += HALF*(W.v.sqr() - TWO*sqr(W.v.x)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
		dFxdU(4,2) -= W.v.x*W.v.y*W.gm1;
		dFxdU(4,3) -= W.v.x*W.v.z*W.gm1;
		dFxdU(4,4) += W.g*W.v.x;
	}

	// y-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fy(void) {
		return (Euler3D_Polytropic_cState(rho*v.y, rho*v.x*v.y, rho*sqr(v.y) + p, rho*v.y*v.z, v.y*H()));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fy(void) const {
		return (Euler3D_Polytropic_cState(rho*v.y, rho*v.x*v.y, rho*sqr(v.y) + p, rho*v.y*v.z, v.y*H()));
	}
	Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho*W.v.y, W.rho*W.v.x*W.v.y, W.rho*sqr(W.v.y) + W.p, W.rho*W.v.y*W.v.z, W.v.y*W.H()));
	}
	void Euler3D_Polytropic_pState::dFydU(DenseMatrix &dFydU) {
		dFydU(0,2) += ONE;
		dFydU(1,0) -= v.x*v.y;
		dFydU(1,1) += v.y;
		dFydU(1,2) += v.x;
		dFydU(2,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.y));
		dFydU(2,1) -= v.x*gm1;
		dFydU(2,2) -= v.y*(g-THREE);
		dFydU(2,3) -= v.z*gm1;
		dFydU(2,4) += gm1;
		dFydU(3,0) -= v.y*v.z;
		dFydU(3,2) += v.z;
		dFydU(3,3) += v.y;
		dFydU(4,0) += HALF*v.y*v.sqr()*(g-TWO) - v.y*p/rho*g*gm1i;
		dFydU(4,1) -= v.x*v.y*gm1;
		dFydU(4,2) += HALF*(v.sqr() - TWO*sqr(v.y)*gm1) + p/rho*g*gm1i;
		dFydU(4,3) -= v.y*v.z*gm1;
		dFydU(4,4) += g*v.y;
	}
	void Euler3D_Polytropic_pState::dFydU(DenseMatrix &dFydU) const {
		dFydU(0,2) += ONE;
		dFydU(1,0) -= v.x*v.y;
		dFydU(1,1) += v.y;
		dFydU(1,2) += v.x;
		dFydU(2,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.y));
		dFydU(2,1) -= v.x*gm1;
		dFydU(2,2) -= v.y*(g-THREE);
		dFydU(2,3) -= v.z*gm1;
		dFydU(2,4) += gm1;
		dFydU(3,0) -= v.y*v.z;
		dFydU(3,2) += v.z;
		dFydU(3,3) += v.y;
		dFydU(4,0) += HALF*v.y*v.sqr()*(g-TWO) - v.y*p/rho*g*gm1i;
		dFydU(4,1) -= v.x*v.y*gm1;
		dFydU(4,2) += HALF*(v.sqr() - TWO*sqr(v.y)*gm1) + p/rho*g*gm1i;
		dFydU(4,3) -= v.y*v.z*gm1;
		dFydU(4,4) += g*v.y;
	}
	void dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_pState &W) {
		dFydU(0,2) += ONE;
		dFydU(1,0) -= W.v.x*W.v.y;
		dFydU(1,1) += W.v.y;
		dFydU(1,2) += W.v.x;
		dFydU(2,0) += HALF*(W.v.sqr()*W.gm1 - TWO*sqr(W.v.y));
		dFydU(2,1) -= W.v.x*W.gm1;
		dFydU(2,2) -= W.v.y*(W.g-THREE);
		dFydU(2,3) -= W.v.z*W.gm1;
		dFydU(2,4) += W.gm1;
		dFydU(3,0) -= W.v.y*W.v.z;
		dFydU(3,2) += W.v.z;
		dFydU(3,3) += W.v.y;
		dFydU(4,0) += HALF*W.v.y*W.v.sqr()*(W.g-TWO) - W.v.y*W.p/W.rho*W.g*W.gm1i;
		dFydU(4,1) -= W.v.x*W.v.y*W.gm1;
		dFydU(4,2) += HALF*(W.v.sqr() - TWO*sqr(W.v.y)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
		dFydU(4,3) -= W.v.y*W.v.z*W.gm1;
		dFydU(4,4) += W.g*W.v.y;
	}

	// z-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) {
		return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z, rho*sqr(v.z) + p, v.z*H()));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) const {
		return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z, rho*sqr(v.z) + p, v.z*H()));
	}
	Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho*W.v.z, W.rho*W.v.x*W.v.z, W.rho*W.v.y*W.v.z, W.rho*sqr(W.v.z) + W.p, W.v.z*W.H()));
	}
	void Euler3D_Polytropic_pState::dFzdU(DenseMatrix &dFzdU) {
		dFzdU(0,3) += ONE;
		dFzdU(1,0) -= v.x*v.z;
		dFzdU(1,1) += v.z;
		dFzdU(1,3) += v.x;
		dFzdU(2,0) -= v.y*v.z;
		dFzdU(2,2) += v.z;
		dFzdU(2,3) += v.y;
		dFzdU(3,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.z));
		dFzdU(3,1) -= v.x*gm1;
		dFzdU(3,2) -= v.y*gm1;
		dFzdU(3,3) -= v.z*(g-THREE);
		dFzdU(3,4) += gm1;
		dFzdU(4,0) += HALF*v.z*v.sqr()*(g-TWO) - v.z*p/rho*g*gm1i;
		dFzdU(4,1) -= v.x*v.z*gm1;
		dFzdU(4,2) -= v.y*v.z*gm1;
		dFzdU(4,3) += HALF*(v.sqr() - TWO*sqr(v.z)*gm1) + p/rho*g*gm1i;
		dFzdU(4,4) += g*v.z;
	}
	void Euler3D_Polytropic_pState::dFzdU(DenseMatrix &dFzdU) const {
		dFzdU(0,3) += ONE;
		dFzdU(1,0) -= v.x*v.z;
		dFzdU(1,1) += v.z;
		dFzdU(1,3) += v.x;
		dFzdU(2,0) -= v.y*v.z;
		dFzdU(2,2) += v.z;
		dFzdU(2,3) += v.y;
		dFzdU(3,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.z));
		dFzdU(3,1) -= v.x*gm1;
		dFzdU(3,2) -= v.y*gm1;
		dFzdU(3,3) -= v.z*(g-THREE);
		dFzdU(3,4) += gm1;
		dFzdU(4,0) += HALF*v.z*v.sqr()*(g-TWO) - v.z*p/rho*g*gm1i;
		dFzdU(4,1) -= v.x*v.z*gm1;
		dFzdU(4,2) -= v.y*v.z*gm1;
		dFzdU(4,3) += HALF*(v.sqr() - TWO*sqr(v.z)*gm1) + p/rho*g*gm1i;
		dFzdU(4,4) += g*v.z;
	}
	void dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_pState &W) {
		dFzdU(0,3) += ONE;
		dFzdU(1,0) -= W.v.x*W.v.z;
		dFzdU(1,1) += W.v.z;
		dFzdU(1,3) += W.v.x;
		dFzdU(2,0) -= W.v.y*W.v.z;
		dFzdU(2,2) += W.v.z;
		dFzdU(2,3) += W.v.y;
		dFzdU(3,0) += HALF*(W.v.sqr()*W.gm1 - TWO*sqr(W.v.z));
		dFzdU(3,1) -= W.v.x*W.gm1;
		dFzdU(3,2) -= W.v.y*W.gm1;
		dFzdU(3,3) -= W.v.z*(W.g-THREE);
		dFzdU(3,4) += W.gm1;
		dFzdU(4,0) += HALF*W.v.z*W.v.sqr()*(W.g-TWO) - W.v.z*W.p/W.rho*W.g*W.gm1i;
		dFzdU(4,1) -= W.v.x*W.v.z*W.gm1;
		dFzdU(4,2) -= W.v.y*W.v.z*W.gm1;
		dFzdU(4,3) += HALF*(W.v.sqr() - TWO*sqr(W.v.z)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
		dFzdU(4,4) += W.g*W.v.z;
	}
	
   /* 
	* Solution variable Jacobians.
	* ----------------------------
	*/

	// dUdW
	void Euler3D_Polytropic_pState::dUdW(DenseMatrix &dUdW) {
		dUdW(0,0) += ONE;
		dUdW(1,0) += v.x;
		dUdW(1,1) += rho;
		dUdW(2,0) += v.y;
		dUdW(2,2) += rho;
		dUdW(3,0) += v.z;
		dUdW(3,3) += rho;
		dUdW(4,0) += HALF*(v.sqr());
		dUdW(4,1) += rho*v.x;
		dUdW(4,2) += rho*v.y;
		dUdW(4,3) += rho*v.z;
		dUdW(4,4) += gm1i;
	}
	void Euler3D_Polytropic_pState::dUdW(DenseMatrix &dUdW) const {
		dUdW(0,0) += ONE;
		dUdW(1,0) += v.x;
		dUdW(1,1) += rho;
		dUdW(2,0) += v.y;
		dUdW(2,2) += rho;
		dUdW(3,0) += v.z;
		dUdW(3,3) += rho;
		dUdW(4,0) += HALF*(v.sqr());
		dUdW(4,1) += rho*v.x;
		dUdW(4,2) += rho*v.y;
		dUdW(4,3) += rho*v.z;
		dUdW(4,4) += gm1i;
	}
	void Euler3D_Polytropic_pState::dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_pState &W) {
		dUdW(0,0) += ONE;
		dUdW(1,0) += W.v.x;
		dUdW(1,1) += W.rho;
		dUdW(2,0) += W.v.y;
		dUdW(2,2) += W.rho;
		dUdW(3,0) += W.v.z;
		dUdW(3,3) += W.rho;
		dUdW(4,0) += HALF*(W.v.sqr());
		dUdW(4,1) += W.rho*W.v.x;
		dUdW(4,2) += W.rho*W.v.y;
		dUdW(4,3) += W.rho*W.v.z;
		dUdW(4,4) += W.gm1i;
	}

	// dWdU
	void Euler3D_Polytropic_pState::dWdU(DenseMatrix &dWdU) {
		dWdU(0,0) += ONE;
		dWdU(1,0) -= v.x/rho;
		dWdU(1,1) += ONE/rho;
		dWdU(2,0) -= v.y/rho;
		dWdU(2,2) += ONE/rho;
		dWdU(3,0) -= v.z/rho;
		dWdU(3,3) += ONE/rho;
		dWdU(4,0) += HALF*(v.sqr())*gm1;
		dWdU(4,1) -= v.x*gm1;
		dWdU(4,2) -= v.y*gm1;
		dWdU(4,3) -= v.z*gm1;
		dWdU(4,4) += gm1;
	}
	void Euler3D_Polytropic_pState::dWdU(DenseMatrix &dWdU) const {
		dWdU(0,0) += ONE;
		dWdU(1,0) -= v.x/rho;
		dWdU(1,1) += ONE/rho;
		dWdU(2,0) -= v.y/rho;
		dWdU(2,2) += ONE/rho;
		dWdU(3,0) -= v.z/rho;
		dWdU(3,3) += ONE/rho;
		dWdU(4,0) += HALF*(v.sqr())*gm1;
		dWdU(4,1) -= v.x*gm1;
		dWdU(4,2) -= v.y*gm1;
		dWdU(4,3) -= v.z*gm1;
		dWdU(4,4) += gm1;
	}
	void Euler3D_Polytropic_pState::dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_pState &W) {
		dWdU(0,0) += ONE;
		dWdU(1,0) -= W.v.x/W.rho;
		dWdU(1,1) += ONE/W.rho;
		dWdU(2,0) -= W.v.y/W.rho;
		dWdU(2,2) += ONE/W.rho;
		dWdU(3,0) -= W.v.z/W.rho;
		dWdU(3,3) += ONE/W.rho;
		dWdU(4,0) += HALF*(W.v.sqr())*W.gm1;
		dWdU(4,1) -= W.v.x*W.gm1;
		dWdU(4,2) -= W.v.y*W.gm1;
		dWdU(4,3) -= W.v.z*W.gm1;
		dWdU(4,4) += W.gm1;
	}



   /* 
	* Eigenvalues
	* ----------- 
	*/

	// x-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_x(void) {
		double c = a();
		return (Euler3D_Polytropic_pState(v.x - c, v.x, v.x, v.x, v.x + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_x(void) const {
		double c = a();
		return (Euler3D_Polytropic_pState(v.x - c, v.x, v.x, v.x, v.x + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_x(const Euler3D_Polytropic_pState &W) {
		double c = W.a();
		return (Euler3D_Polytropic_pState(W.v.x - c, W.v.x, W.v.x, W.v.x, W.v.x + c));
	}
	double Euler3D_Polytropic_pState::lambda_x(int index) {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.x-a());
			case 2 :
				return (v.x);
			case 3 :
				return (v.x);
			case 4 :
				return (v.x);
			case 5 :
				return (v.x+a());
			default:
				return (v.x);
		}
	}
	double Euler3D_Polytropic_pState::lambda_x(int index) const {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.x-a());
			case 2 :
				return (v.x);
			case 3 :
				return (v.x);
			case 4 :
				return (v.x);
			case 5 :
				return (v.x+a());
			default:
				return (v.x);
		}
	}

	// y-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_y(void) {
		double c = a();
		return (Euler3D_Polytropic_pState(v.y - c, v.y, v.y, v.y, v.y + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_y(void) const {
		double c = a();
		return (Euler3D_Polytropic_pState(v.y - c, v.y, v.y, v.y, v.y + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_y(const Euler3D_Polytropic_pState &W) {
		double c = W.a();
		return (Euler3D_Polytropic_pState(W.v.y - c, W.v.y, W.v.y, W.v.y, W.v.y + c));
	}
	double Euler3D_Polytropic_pState::lambda_y(int index) {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.y-a());
			case 2 :
				return (v.y);
			case 3 :
				return (v.y);
			case 4 :
				return (v.y);
			case 5 :
				return (v.y+a());
			default:
				return (v.y);
		};
	}
	double Euler3D_Polytropic_pState::lambda_y(int index) const {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.y-a());
			case 2 :
				return (v.y);
			case 3 :
				return (v.y);
			case 4 :
				return (v.y);
			case 5 :
				return (v.y+a());
			default:
				return (v.y);
		}
	}
	
	
	// z-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_z(void) {
		double c = a();
		return (Euler3D_Polytropic_pState(v.z - c, v.z, v.z, v.z, v.z + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_z(void) const {
		double c = a();
		return (Euler3D_Polytropic_pState(v.z - c, v.z, v.z, v.z, v.z + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda_z(const Euler3D_Polytropic_pState &W) {
		double c = W.a();
		return (Euler3D_Polytropic_pState(W.v.z - c, W.v.z, W.v.z, W.v.z, W.v.z + c));
	}
	double Euler3D_Polytropic_pState::lambda_z(int index) {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.z-a());
			case 2 :
				return (v.z);
			case 3 :
				return (v.z);
			case 4 :
				return (v.z);
			case 5 :
				return (v.z+a());
			default:
				return (v.z);
		}
	}
	double Euler3D_Polytropic_pState::lambda_z(int index) const {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (v.z-a());
			case 2 :
				return (v.z);
			case 3 :
				return (v.z);
			case 4 :
				return (v.z);
			case 5 :
				return (v.z+a());
			default:
				return (v.z);
		}
	}



   /*
	* Conserved right eigenvectors
	* ----------------------------
	*/

	// x-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_x(const int &index){
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a()));
			}		
		}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_x(const int &index) const{
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a()));
		}		
	}

	// y-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_y(const int &index) {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y-a(), v.z, h()-v.y*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr()));
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, rho, ZERO, ZERO, rho*v.x));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y+a(), v.z, h()+v.y*a()));
		}		
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_y(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y-a(), v.z, h()-v.y*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr()));
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, rho, ZERO, ZERO, rho*v.x));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y+a(), v.z, h()+v.y*a()));
		}		
	}

	// z-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_z(const int &index) {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z-a(), h()-v.z*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr()));
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, rho, ZERO, ZERO, rho*v.x));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z+a(), h()+v.z*a()));
		}		
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_z(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z-a(), h()-v.z*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr()));
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, rho, ZERO, ZERO, rho*v.x));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z+a(), h()+v.z*a()));
		}		
	}

   /*
	* Primitive left eigenvectors
	* ---------------------------
	*/

	// x-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_x(const int &index) {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, -rho/(TWO*a()), ZERO, ZERO, HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, rho/(TWO*a()), ZERO, ZERO, HALF/a2()));
		}		
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_x(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, -rho/(TWO*a()), ZERO, ZERO, HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, rho/(TWO*a()), ZERO, ZERO, HALF/a2()));
		}		
	}

	// y-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_y(const int &index) {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, -rho/(TWO*a()), ZERO, HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, rho/(TWO*a()), ZERO, HALF/a2()));
		}		
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_y(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, -rho/(TWO*a()), ZERO, HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, rho/(TWO*a()), ZERO, HALF/a2()));
		}		
	}

	// z-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_z(const int &index) {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, -rho/(TWO*a()), HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, rho/(TWO*a()), HALF/a2()));
			}		
		}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_z(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, -rho/(TWO*a()), HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, rho/(TWO*a()), HALF/a2()));	
		}
					
	}

   /*
	* Primitive right eigenvector
	* ---------------------------
	*/

	// x-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_x(const int &index) {
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, -a()/rho, ZERO, ZERO, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, a()/rho, ZERO, ZERO, sqr(a())));
		}		
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_x(const int &index) const {
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, -a()/rho, ZERO, ZERO, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, a()/rho, ZERO, ZERO, sqr(a())));
		}		
	}

	// y-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_y(const int &index) {
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, -a()/rho, ZERO, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, a()/rho, ZERO, sqr(a())));
		}		
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_y(const int &index) const {
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, -a()/rho, ZERO, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ZERO, ONE, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, a()/rho, ZERO, sqr(a())));
		}		
	}


	// z-direction
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_z(const int &index) {
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, -a()/rho, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, a()/rho, sqr(a())));
		}		
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp_z(const int &index) const{
		switch(index) {
			case 1: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, -a()/rho, sqr(a())));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, ZERO)); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
			case 4: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, ONE, ZERO, ZERO));
			case 5: 
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, a()/rho, sqr(a())));
		}		
	}


/********************************************************
* Routine: Reflect                                     *
*                                                      *
* This function returns the reflected solution state   *
* in a given direction given the primitive solution    *
* variables and the unit normal vector in the          *
* direction of interest.                               *
*                                                      *
********************************************************/
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::Reflect(const Euler3D_Polytropic_pState &W, const Vector3D &norm_dir) {
	
	Vector3D ur_norm, ur_tang, vr_tot;
	
	Euler3D_Polytropic_pState Temp;
	Temp.Copy(W);
	
	ur_norm = dot(W.v, norm_dir)*norm_dir;
	ur_tang = W.v - ur_norm;
	
	ur_norm = -ur_norm;
	vr_tot = ur_norm + ur_tang;
	
	Temp.v = vr_tot;
	
    return (Temp);
}

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::Moving_Wall(
									 const Euler3D_Polytropic_pState &Win,
									 const Euler3D_Polytropic_pState &Wout,
									 const Vector3D &norm_dir, 
									 const Vector3D &wall_velocity,
									 const Vector3D &pressure_gradient,
									 const int &TEMPERATURE_BC_FLAG) {
	
	Euler3D_Polytropic_pState Temp;
	Temp.Copy(Win);
	
	if(wall_velocity == Vector3D_ZERO){
		Temp.v = -Win.v;
	}else{
		
		double  Wall_velocity_tang ;
		Vector3D ur_norm, ur_tang, vr_tot, uw_tang;
		
		
		ur_norm = dot(Win.v, norm_dir)*norm_dir;
		ur_tang = Win.v - ur_norm;
		
		
		uw_tang = wall_velocity - dot(norm_dir,wall_velocity)*norm_dir;
		
		ur_norm = -ur_norm;
		ur_tang = 2.0*uw_tang - ur_tang;
		vr_tot = ur_norm + ur_tang;
		
		Temp.v = vr_tot;
	}
	
	//  Fixed Wall Temperature or constant extrapolation for Adiabatic
	if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
		if (pressure_gradient != Vector3D_ZERO){
			Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
		}else{
			
			Temp.rho = Win.p/(Temp.Rtot()*Wout.T());
		}
		
	}
	
	return (Temp);
	
}

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::No_Slip(
										 const Euler3D_Polytropic_pState &Win,
										 const Euler3D_Polytropic_pState &Wout,
										 const Vector3D &norm_dir,
										 const Vector3D &pressure_gradient,
										 const int &TEMPERATURE_BC_FLAG) {
	
	return (Moving_Wall(Win, Wout, norm_dir,Vector3D_ZERO,pressure_gradient,TEMPERATURE_BC_FLAG));
	
}







/* ------------------------------------------------------------------------- *
 *				 Euler3D_Polytropic_cState subroutines						 *
 * ------------------------------------------------------------------------- */

   /*
	* Constructors
	* ------------
	*/

	// Creation constructor
	Euler3D_Polytropic_cState::Euler3D_Polytropic_cState() {
		rho = DENSITY_STDATM;	rhov.zero();	E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
	}

	// Copy constructor
	Euler3D_Polytropic_cState::Euler3D_Polytropic_cState(const Euler3D_Polytropic_cState &U) {
		rho = U.rho;	rhov = U.rhov;		E = U.E;
	}

	// Assignment constructors
	Euler3D_Polytropic_cState::Euler3D_Polytropic_cState(const double &d, 
							  const Vector3D &dv, 
							  const double &Etotal) {
		rho = d;	rhov = dv;		E = Etotal;
	}

	Euler3D_Polytropic_cState::Euler3D_Polytropic_cState(const double &d, 
							  const double &dvx, const double &dvy, const double &dvz, 
							  const double &Etotal) {
		rho = d;	
		rhov.x = dvx;	rhov.y = dvy;	rhov.z = dvz;
		E = Etotal;
	}

	Euler3D_Polytropic_cState::Euler3D_Polytropic_cState(const Euler3D_Polytropic_pState &W) {
		rho = W.rho;	rhov = W.rhov();	E = W.E();
	}





   /*
	* Set static variables 
	* -------------------- 
	*/

	int Euler3D_Polytropic_cState::num_vars = NUM_VAR_EULER3D;
	double Euler3D_Polytropic_cState::g = GAMMA_AIR;
	double Euler3D_Polytropic_cState::gm1 = GAMMA_AIR-ONE;
	double Euler3D_Polytropic_cState::gm1i = ONE/(GAMMA_AIR-ONE);
	double Euler3D_Polytropic_cState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);

	// set gas constants.
	void Euler3D_Polytropic_cState::setgas(void) {
		g = GAMMA_AIR;
		R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		gm1 = g - ONE;
		gm1i = ONE/gm1;
	}
	void Euler3D_Polytropic_cState::setgas(char *string_ptr) {
		if (strcmp(string_ptr, "AIR") == 0) {
			g = GAMMA_AIR;
			R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		} else if (strcmp(string_ptr, "A") == 0) {
			g = GAMMA_A;
			R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
		} else if (strcmp(string_ptr, "CO") == 0) {
			g = GAMMA_CO;
			R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
		} else if (strcmp(string_ptr, "CO2") == 0) {
			g = GAMMA_CO2;
			R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
		} else if (strcmp(string_ptr, "CH4") == 0) {
			g = GAMMA_CH4;
			R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
		} else if (strcmp(string_ptr, "H") == 0) {
			g = GAMMA_H;
			R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
		} else if (strcmp(string_ptr, "H2") == 0) {
			g = GAMMA_H2;
			R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
		} else if (strcmp(string_ptr, "HE") == 0) {
			g = GAMMA_HE;
			R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
		} else if (strcmp(string_ptr, "H2O") == 0) {
			g = GAMMA_H2O;
			R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
		} else if (strcmp(string_ptr, "N2") == 0) {
			g = GAMMA_N2;
			R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
		} else if (strcmp(string_ptr, "O") == 0) {
			g = GAMMA_O;
			R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
		} else if (strcmp(string_ptr, "O2") == 0) {
			g = GAMMA_O2;
			R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
		} else if (strcmp(string_ptr, "e") == 0) {
			g = GAMMA_e;
			R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
		} else {
			g = GAMMA_AIR;
			R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
		}
		gm1 = g - ONE;
	    gm1i = ONE/gm1;
	}


// Flow velocity. 
Vector3D Euler3D_Polytropic_cState::v(void) {
    return (rhov/rho);
}
Vector3D Euler3D_Polytropic_cState::v(void) const {
    return (rhov/rho);
}
double Euler3D_Polytropic_cState::v(const Vector3D &n) {
    return ((rhov*n)/rho);
}
double Euler3D_Polytropic_cState::v(const Vector3D &n) const {
    return ((rhov*n)/rho);
}

// Pressure. 
double Euler3D_Polytropic_cState::p(void) {
    return (gm1*(E - HALF*rhov.sqr()/rho));
}
double Euler3D_Polytropic_cState::p(void) const {
    return (gm1*(E - HALF*rhov.sqr()/rho));
}

// Temperature. 
double Euler3D_Polytropic_cState::T(void) {
    return (gm1*(E - HALF*rhov.sqr()/rho)/(rho*R));
}
double Euler3D_Polytropic_cState::T(void) const {
    return (gm1*(E - HALF*rhov.sqr()/rho)/(rho*R));
}

// Specific internal energy. 
double Euler3D_Polytropic_cState::e(void) {
    return (E/rho - HALF*rhov.sqr()/sqr(rho));
}
double Euler3D_Polytropic_cState::e(void) const {
    return (E/rho - HALF*rhov.sqr()/sqr(rho));
}

// Specific enthalpy. 
double Euler3D_Polytropic_cState::h(void) {
    return (g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho));
}
double Euler3D_Polytropic_cState::h(void) const {
    return (g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho));
}

// Total enthalpy. 
double Euler3D_Polytropic_cState::H(void) {
	return (g*E - gm1*HALF*rhov.sqr()/rho);
}
double Euler3D_Polytropic_cState::H(void) const {
	return (g*E - gm1*HALF*rhov.sqr()/rho);
}

// Sound speed. 
double Euler3D_Polytropic_cState::a(void) {
    return (sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))));
}
double Euler3D_Polytropic_cState::a(void) const {
    return (sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))));
}

// Sound speed squared. 
double Euler3D_Polytropic_cState::a2(void) {
    return (g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)));
}
double Euler3D_Polytropic_cState::a2(void) const {
    return (g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)));
}

// Mach number. 
double Euler3D_Polytropic_cState::M(void) {
    return (abs(rhov)/(rho*sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}
double Euler3D_Polytropic_cState::M(void) const {
    return (abs(rhov)/(rho*sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

// Specific entropy. 
double Euler3D_Polytropic_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*rhov.sqr()/rho)/pow(rho, g)));
}
double Euler3D_Polytropic_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*rhov.sqr()/rho)/pow(rho, g)));
}

// Stagnation temperature. 
double Euler3D_Polytropic_cState::To(void) {
    return ((gm1*(E - HALF*rhov.sqr()/rho)/(rho*R))*
			(ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}
double Euler3D_Polytropic_cState::To(void) const {
    return ((gm1*(E - HALF*rhov.sqr()/rho)/(rho*R))*
			(ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

// Stagnation pressure. 
double Euler3D_Polytropic_cState::po(void) {
    return ((gm1*(E - HALF*rhov.sqr()/rho))*
            pow(ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*
                                       (E/rho -HALF*rhov.sqr()/sqr(rho))), g*gm1i));
}
double Euler3D_Polytropic_cState::po(void) const {
    return ((gm1*(E - HALF*rhov.sqr()/rho))*
            pow(ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*
                                       (E/rho - HALF*rhov.sqr()/sqr(rho))), g*gm1i));
}

// Stagnation sound speed. 
double Euler3D_Polytropic_cState::ao(void) {
	return (sqrt((g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))*
				 (ONE+HALF*gm1*rhov.sqr()/
				  (rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))))));
}
double Euler3D_Polytropic_cState::ao(void) const {
	return (sqrt((g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))*
				 (ONE+HALF*gm1*rhov.sqr()/
				  (rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))))));
}


// Stagnation enthalpy. 
double Euler3D_Polytropic_cState::ho(void) {
	return ((g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho))*
			(ONE+HALF*gm1*rhov.sqr()/
			 (rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}
double Euler3D_Polytropic_cState::ho(void) const {
    return ((g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho))*
            (ONE+HALF*gm1*rhov.sqr()/
             (rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}


   /*
	* Primitive solution state. 
	* -------------------------
	*/
	Euler3D_Polytropic_pState Euler3D_Polytropic_cState::W(void) {
		return (Euler3D_Polytropic_pState(rho, v(), p()));
	}

	Euler3D_Polytropic_pState Euler3D_Polytropic_cState::W(void) const {
		return (Euler3D_Polytropic_pState(rho, v(), p()));
	}

	Euler3D_Polytropic_pState Euler3D_Polytropic_cState::W(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_pState(U.rho, U.v(), U.p()));
	}


   /* 
	* Fluxes and Jacobians
	* --------------------
	*/

	// x-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(void) {
		return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(void) const {
		return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rhov.x, sqr(U.rhov.x)/U.rho + U.p(), U.rhov.x*U.rhov.y/U.rho, U.rhov.x*U.rhov.z/U.rho, U.rhov.x*U.H()/U.rho));
	}
	void Euler3D_Polytropic_cState::dFxdU(DenseMatrix &dFxdU) {
		W().dFxdU(dFxdU);
	}
	void Euler3D_Polytropic_cState::dFxdU(DenseMatrix &dFxdU) const {
		W().dFxdU(dFxdU);
	}
	void Euler3D_Polytropic_cState::dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_cState &U) {
		U.W().dFxdU(dFxdU);
	}

	// y-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(void) {
		return (Euler3D_Polytropic_cState(rhov.y, rhov.x*rhov.y/rho, sqr(rhov.y)/rho + p(), rhov.y*rhov.z/rho, rhov.y*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(void) const {
		return (Euler3D_Polytropic_cState(rhov.y, rhov.x*rhov.y/rho, sqr(rhov.y)/rho + p(), rhov.y*rhov.z/rho, rhov.y*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rhov.y, U.rhov.x*U.rhov.y/U.rho, sqr(U.rhov.y)/U.rho + U.p(), U.rhov.y*U.rhov.z/U.rho, U.rhov.y*U.H()/U.rho));
	}
	void Euler3D_Polytropic_cState::dFydU(DenseMatrix &dFydU) {
		W().dFydU(dFydU);
	}
	void Euler3D_Polytropic_cState::dFydU(DenseMatrix &dFydU) const {
		W().dFydU(dFydU);
	}
	void Euler3D_Polytropic_cState::dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_cState &U) {
		U.W().dFydU(dFydU);
	}

	// z-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(void) {
		return (Euler3D_Polytropic_cState(rhov.z, rhov.x*rhov.z/rho, rhov.y*rhov.z/rho, sqr(rhov.z)/rho + p(), rhov.z*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(void) const {
		return (Euler3D_Polytropic_cState(rhov.z, rhov.x*rhov.z/rho, rhov.y*rhov.z/rho, sqr(rhov.z)/rho + p(), rhov.z*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rhov.z, U.rhov.x*U.rhov.z/U.rho, U.rhov.y*U.rhov.z/U.rho, sqr(U.rhov.z)/U.rho + U.p(), U.rhov.z*U.H()/U.rho));
	}
	void Euler3D_Polytropic_cState::dFzdU(DenseMatrix &dFzdU) {
		W().dFzdU(dFzdU);
	}
	void Euler3D_Polytropic_cState::dFzdU(DenseMatrix &dFzdU) const {
		W().dFzdU(dFzdU);
	}
	void Euler3D_Polytropic_cState::dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_cState &U) {
		U.W().dFzdU(dFzdU);
	}

   /*
	* Solution variable Jacobian.
	* ---------------------------
	*/

	// dUdW
	void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW) {
		W().dUdW(dUdW);
	}
	void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW) const {
		W().dUdW(dUdW);
	}
	void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_cState &U) {
		U.W().dUdW(dUdW);
	}

	// dWdU
	void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU) {
		W().dWdU(dWdU);
	}
	void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU) const {
		W().dWdU(dWdU);
	}
	void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_cState &U) {
		U.W().dWdU(dWdU);
	}

   /*
	* Operators.
	* ----------
	*/

	// Index operator. 
	double& Euler3D_Polytropic_cState::operator[](int index) {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (rho);
			case 2 :
				return (rhov.x);
			case 3 :
				return (rhov.y);
			case 4 :
				return (rhov.z);
			case 5 :
				return (E);
			default:
				return (rho);
		}	
	}
	const double& Euler3D_Polytropic_cState::operator[](int index) const {
		assert( index >= 1 && index <= NUM_VAR_EULER3D );
		switch(index) {
			case 1 :
				return (rho);
			case 2 :
				return (rhov.x);
			case 3 :
				return (rhov.y);
			case 4 :
				return (rhov.z);
			case 5 :
				return (E);
			default:
				return (rho);
		}	
	}

	// Binary arithmetic operators.
	Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho+U2.rho,U1.rhov+U2.rhov,U1.E+U2.E));
	}

	Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho-U2.rho,U1.rhov-U2.rhov,U1.E-U2.E));
	}

	double operator *(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (U1.rho*U2.rho+U1.rhov*U2.rhov+U1.E*U2.E);
	}

	Euler3D_Polytropic_cState operator *(const Euler3D_Polytropic_cState &U, const double &a) {
		return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
	}

	Euler3D_Polytropic_cState operator *(const double &a, const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
	}
	Euler3D_Polytropic_cState operator /(const Euler3D_Polytropic_cState &U, const double &a) {
		return (Euler3D_Polytropic_cState(U.rho/a,U.rhov/a,U.E/a));
	}

	Euler3D_Polytropic_cState operator ^(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho*U2.rho,
										  U1.rhov.x*U2.rhov.x,
										  U1.rhov.y*U2.rhov.y,
										  U1.rhov.z*U2.rhov.z,
										  U1.E*U2.E));
	}

	// Unary arithmetic operators.
	Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rho,U.rhov,U.E));
	}
	Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(-U.rho,-U.rhov,-U.E));
	}

	// Shortcut arithmetic operators.
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator +=(const Euler3D_Polytropic_cState &U) {
		rho += U.rho;
		rhov += U.rhov;
		E += U.E;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator -=(const Euler3D_Polytropic_cState &U) {
		rho -= U.rho;
		rhov -= U.rhov;
		E -= U.E;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator *=(const double &a) {
		rho *= a;
		rhov.x *= a;		rhov.y *= a;		rhov.z *= a;
		E *= a;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator /=(const double &a) {
		rho /= a;
		rhov.x /= a;		rhov.y /= a;		rhov.z /= a;
		E /= a;
		return *this;
	}


	// Relational operators.
	int operator ==(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E);
	}

	int operator !=(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E);
	}


	// Input-output operators.
	ostream& operator << (ostream &out_file, const Euler3D_Polytropic_cState &U) {
		out_file.setf(ios::scientific);
		out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " "
			<< U.rhov.z << " " << U.E;
		out_file.unsetf(ios::scientific);
		return (out_file);
	}
	istream& operator >> (istream &in_file, Euler3D_Polytropic_cState &U) {
		in_file.setf(ios::skipws);
		in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
		in_file.unsetf(ios::skipws);
		return (in_file);
	}





/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::RoeAverage(const Euler3D_Polytropic_pState &Wl, const Euler3D_Polytropic_pState &Wr) {

    double hl, hr, sqrt_rhol, sqrt_rhor;
    double da, ua, va,wa, pa, aa2, ha, ga, gam1;

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */

    hl = Wl.h();
    hr = Wr.h();
    sqrt_rhol = sqrt(Wl.rho);
    sqrt_rhor = sqrt(Wr.rho);

    /* Determine the appropriate Roe averages. */
    ga = Wl.g;
    gam1 = Wl.gm1;
    da = sqrt_rhol*sqrt_rhor;
    ua = (sqrt_rhol*Wl.v.x+sqrt_rhor*Wr.v.x)/(sqrt_rhol+sqrt_rhor);
    va = (sqrt_rhol*Wl.v.y+sqrt_rhor*Wr.v.y)/(sqrt_rhol+sqrt_rhor);
    wa = (sqrt_rhol*Wl.v.z+sqrt_rhor*Wr.v.z)/(sqrt_rhol+sqrt_rhor);

    ha = (sqrt_rhol*hl+sqrt_rhor*hr)/(sqrt_rhol+sqrt_rhor);
    aa2 = gam1*(ha-HALF*(sqr(ua)+sqr(va)+sqr(wa)));
    pa = da*aa2/ga;

    /* Return the Roe-averged state. */

    return (Euler3D_Polytropic_pState(da, ua, va, wa, pa));
}

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::RoeAverage(const Euler3D_Polytropic_cState &Ul, const Euler3D_Polytropic_cState &Ur) {
	return (RoeAverage(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_x (Harten-Lax-van Leer flux         *
 *                      function, x-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_x(const Euler3D_Polytropic_pState &Wl,
	      	          const Euler3D_Polytropic_pState &Wr) {

    double wavespeed_l, wavespeed_r;
   
    Euler3D_Polytropic_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Euler3D_Polytropic_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER3D],
                      lambdas_a[NUM_VAR_EULER3D]);
 
    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fx();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fx();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
	      	          const Euler3D_Polytropic_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
}


Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_y(const Euler3D_Polytropic_pState &Wl,
									 const Euler3D_Polytropic_pState &Wr) {
	
    double wavespeed_l, wavespeed_r;
	
    Euler3D_Polytropic_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Euler3D_Polytropic_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);
	
	/* Evaluate the jumps in the conserved solution states. */
	
    dUrl = Wr.U()-Wl.U();
	
    /* Evaluate the left, right, and average state eigenvalues. */
	
    lambdas_l = Wl.lambda_y();
    lambdas_r = Wr.lambda_y();
    lambdas_a = Wa.lambda_y();
	
    /* Determine the intermediate state flux. */
	
    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER3D],
                      lambdas_a[NUM_VAR_EULER3D]);
	
    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fy();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fy();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fy()-wavespeed_l*Wr.Fy())
				   +(wavespeed_l*wavespeed_r)*dUrl)/
		(wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);
}
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_y(const Euler3D_Polytropic_cState &Ul,
									 const Euler3D_Polytropic_cState &Ur) {
	return (FluxHLLE_y(Ul.W(), Ur.W()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_z(const Euler3D_Polytropic_pState &Wl,
									 const Euler3D_Polytropic_pState &Wr) {
	
    double wavespeed_l, wavespeed_r;
	
    Euler3D_Polytropic_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Euler3D_Polytropic_cState Flux, dUrl;
    
    /* Evaluate the Roe-average primitive solution state. */
    
    Wa = RoeAverage(Wl, Wr);
	
	/* Evaluate the jumps in the conserved solution states. */
	
    dUrl = Wr.U()-Wl.U();
	
    /* Evaluate the left, right, and average state eigenvalues. */
	
    lambdas_l = Wl.lambda_z();
    lambdas_r = Wr.lambda_z();
    lambdas_a = Wa.lambda_z();
	
    /* Determine the intermediate state flux. */
	
    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed_r = max(lambdas_r[NUM_VAR_EULER3D],
                      lambdas_a[NUM_VAR_EULER3D]);
	
    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fz();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fz();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fz()-wavespeed_l*Wr.Fz())
				   +(wavespeed_l*wavespeed_r)*dUrl)/
		(wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_z(const Euler3D_Polytropic_cState &Ul,
									 const Euler3D_Polytropic_cState &Ur) {
	return (FluxHLLE_z(Ul.W(), Ur.W()));
}

/*********************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux         *
 *                      function, n-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the so-called         *
 * Harten-Lax-van Leer approximation to specify the      *
 * intermediate state fluxes in terms of the rotated     *
 * solution states.  See Harten, Lax, van Leer (1983).   *
 *                                                       *
 *********************************************************/
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_n(const Euler3D_Polytropic_pState &Wl,
	      	          const Euler3D_Polytropic_pState &Wr,
                          const Vector3D &norm_dir) {

	double Wl_ur_norm, Wl_ur_tang;
	double Wr_ur_norm, Wr_ur_tang ;
	double Wr_ur_tang_z;
	
	Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
	Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
	Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
	Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
	
	//solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
	Euler3D_Polytropic_pState Wl_rotated, Wr_rotated;
	Euler3D_Polytropic_cState Flux, Flux_rotated;
	
	
	/* Apply the frame rotation and evaluate left and right
		solution states in the local rotated frame defined
		by the unit normal vector. */
	Wl_rotated.Copy(Wl);
	Wr_rotated.Copy(Wr);
	
	// Left state velocity in rotated frame
	Wl_ur_norm = dot(Wl.v, norm_dir);
	Wl_ur_tang = abs(Wl.v - Wl_ur_norm*norm_dir);
	Wl_ur_tang_vector = (Wl.v - Wl_ur_norm*norm_dir);
	if(Wl_ur_tang != ZERO){
		Wl_ur_tang_unit_vector =  Wl_ur_tang_vector/Wl_ur_tang;
	}else{
		Wl_ur_tang_unit_vector= Vector3D_ZERO;
	}
	
	Wl_rotated.rho = Wl.rho;
	Wl_rotated.v.x = Wl_ur_norm ;
	Wl_rotated.v.y = Wl_ur_tang;
	Wl_rotated.v.z = ZERO;
	Wl_rotated.p = Wl.p;
	// Right state velocity in rotated frame
	Wr_ur_norm = dot(Wr.v, norm_dir);
	Wr_ur_tang_vector = Wr.v - Wr_ur_norm*norm_dir;
	Wr_ur_tang = abs(Wr.v - Wr_ur_norm*norm_dir);
	if( Wr_ur_tang != ZERO){
		Wr_ur_tang_unit_vector =  Wr_ur_tang_vector/Wr_ur_tang ;
	}else{
		Wr_ur_tang_unit_vector= Vector3D_ZERO;  
	}
	
	Wr_rotated.rho = Wr.rho;
	Wr_rotated.v.x = Wr_ur_norm;
	Wr_rotated.v.y = dot( Wr_ur_tang_vector, Wl_ur_tang_unit_vector);
	Wr_rotated.v.z = abs( Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector);
	Wr_rotated.p = Wr.p;
	
	
	Wr_ur_tang_z = abs(Wr_ur_tang_vector-Wr_rotated.v.y* Wl_ur_tang_unit_vector);
	Wr_ur_tang_z_vector = Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector;
	if(Wr_ur_tang_z !=ZERO){
		Wr_ur_tang_z_unit_vector = Wr_ur_tang_z_vector /Wr_ur_tang_z ;
	}else{
		Wr_ur_tang_z_unit_vector = Vector3D_ZERO;
		
	}
	
	
    /* Evaluate the intermediate state solution 
		flux in the rotated frame. */
	
	Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);
	
	/* Rotate back to the original Cartesian reference
		frame and return the solution flux. */
	
	Flux.Copy(Flux_rotated);
	
	
	
	Flux_rotated_x = Flux.rhov.x*norm_dir;
	Flux_rotated_tang_y = Flux.rhov.y* Wl_ur_tang_unit_vector ;
	Flux_rotated_tang_z = Flux.rhov.z* Wr_ur_tang_z_unit_vector;
	
	
	Flux.rhov =  Flux_rotated_x + Flux_rotated_tang_y+ Flux_rotated_tang_z;
	
	//Flux.zero_non_sol();
	
	
	return (Flux);
	
}


Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
	      	          const Euler3D_Polytropic_cState &Ur,
                          const Vector3D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}


/********************************************************
* Routine: HartenFixPos (Harten Entropy Fix)           *
*                                                      *
* This function returns the positive parts of the      *
* corrected elemental wave speeds or eigenvalues       *
* according to the entropy fix of Harten (1983).       *
*                                                       *
********************************************************/
Euler3D_Polytropic_pState HartenFixPos(
											 const Euler3D_Polytropic_pState &lambdas_a,
											 const Euler3D_Polytropic_pState &lambdas_l,
											 const Euler3D_Polytropic_pState &lambdas_r) {
	
	Euler3D_Polytropic_pState NEW;
	NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
	NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
	NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
	NEW.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
	
	NEW.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);	
	return (NEW);
}


/********************************************************
* Routine: HartenFixNeg (Harten Entropy Fix)           *
*                                                      *
* This function returns the negative parts of the      *
* corrected elemental wave speeds or eigenvalues       *
* according to the entropy fix of Harten (1983).       *
*                                                      *
********************************************************/
Euler3D_Polytropic_pState HartenFixNeg(
											 const Euler3D_Polytropic_pState &lambdas_a,
											 const Euler3D_Polytropic_pState &lambdas_l,
											 const Euler3D_Polytropic_pState &lambdas_r) {
	
	Euler3D_Polytropic_pState NEW;
	
	NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
	NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
	NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
	NEW.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
	NEW.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
		
	return (NEW);
}


// Flux Roe -- based on Harten fix 
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_x(
																		   const  Euler3D_Polytropic_pState &Wl,  
																		   const  Euler3D_Polytropic_pState &Wr){
	
	Euler3D_Polytropic_pState Wa, dWrl, wavespeeds, 
	lambdas_l, lambdas_r, lambdas_a;
	Euler3D_Polytropic_cState Flux;
	
	
	/* Evaluate the Roe-average primitive solution state. */    
	Wa = RoeAverage(Wl, Wr);
	
	/* Evaluate the jumps in the primitive solution states. */
	dWrl = Wr-Wl;
	
	/* Evaluate the left, right, and average state eigenvalues. */
	lambdas_l = Wl.lambda_x();
	lambdas_r = Wr.lambda_x();
	lambdas_a = Wa.lambda_x();
	
	/* Determine the intermediate state flux. */
	if (Wa.v.x >= ZERO) {
		Flux = Wl.Fx();   
		
		wavespeeds = HartenFixNeg(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		
		
		
		for (int i=1 ; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] < ZERO) {
				Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
				
				
			}
		} 
	} else {
		
		Flux = Wr.Fx();
		
		wavespeeds = HartenFixPos(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		for (int i=1; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] > ZERO) {
				Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
				
				
			}
		} 
	} 
    
	/* Return solution flux. */    
	return (Flux);    
	
}
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_x(const  Euler3D_Polytropic_cState &Ul,  
															   const  Euler3D_Polytropic_cState &Ur) {
	return (FluxRoe_x(Ul.W(),Ur.W()));
}

// Flux Roe -- based on Harten fix 
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_y(
															   const  Euler3D_Polytropic_pState &Wl,  
															   const  Euler3D_Polytropic_pState &Wr){
	
	Euler3D_Polytropic_pState Wa, dWrl, wavespeeds, 
	lambdas_l, lambdas_r, lambdas_a;
	Euler3D_Polytropic_cState Flux;
	
	
	/* Evaluate the Roe-average primitive solution state. */    
	Wa = RoeAverage(Wl, Wr);
	
	/* Evaluate the jumps in the primitive solution states. */
	dWrl = Wr-Wl;
	
	/* Evaluate the left, right, and average state eigenvalues. */
	lambdas_l = Wl.lambda_y();
	lambdas_r = Wr.lambda_y();
	lambdas_a = Wa.lambda_y();
	
	/* Determine the intermediate state flux. */
	if (Wa.v.x >= ZERO) {
		Flux = Wl.Fy();   
		
		wavespeeds = HartenFixNeg(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		
		
		
		for (int i=1 ; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] < ZERO) {
				Flux += wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
				
				
			}
		} 
	} else {
		
		Flux = Wr.Fy();
		
		wavespeeds = HartenFixPos(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		for (int i=1; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] > ZERO) {
				Flux -= wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i); 
				
				
			}
		} 
	} 
    
	/* Return solution flux. */    
	return (Flux);    
}
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_y(const  Euler3D_Polytropic_cState &Ul,  
															   const  Euler3D_Polytropic_cState &Ur) {
	return (FluxRoe_y(Ul.W(),Ur.W()));
}

// Flux Roe -- based on Harten fix 
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_z(const  Euler3D_Polytropic_pState &Wl,  
															   const  Euler3D_Polytropic_pState &Wr){
	
	Euler3D_Polytropic_pState Wa, dWrl, wavespeeds, 
	lambdas_l, lambdas_r, lambdas_a;
	Euler3D_Polytropic_cState Flux;
	
	
	/* Evaluate the Roe-average primitive solution state. */    
	Wa = RoeAverage(Wl, Wr);
	
	/* Evaluate the jumps in the primitive solution states. */
	dWrl = Wr-Wl;
	
	/* Evaluate the left, right, and average state eigenvalues. */
	lambdas_l = Wl.lambda_z();
	lambdas_r = Wr.lambda_z();
	lambdas_a = Wa.lambda_z();
	
	/* Determine the intermediate state flux. */
	if (Wa.v.x >= ZERO) {
		Flux = Wl.Fz();   
		
		wavespeeds = HartenFixNeg(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		
		
		
		for (int i=1 ; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] < ZERO) {
				Flux += wavespeeds[i]*(Wa.lp_z(i)*dWrl)*Wa.rc_z(i);
				
				
			}
		} 
	} else {
		
		Flux = Wr.Fz();
		
		wavespeeds = HartenFixPos(lambdas_a,
								  lambdas_l,
								  lambdas_r);
		
		for (int i=1; i <= Wl.num_vars; i++) {
			if (wavespeeds[i] > ZERO) {
				Flux -= wavespeeds[i]*(Wa.lp_z(i)*dWrl)*Wa.rc_z(i); 
				
				
			}
		} 
	} 
    
	/* Return solution flux. */    
	return (Flux);    
    
	
}
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_z(const  Euler3D_Polytropic_cState &Ul,  
															   const  Euler3D_Polytropic_cState &Ur) {
	return (FluxRoe_y(Ul.W(),Ur.W()));
}


Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_n(const Euler3D_Polytropic_pState &Wl,
															   const Euler3D_Polytropic_pState &Wr,
															   const Vector3D &norm_dir){
	
	double Wl_ur_norm, Wl_ur_tang;
	double Wr_ur_norm, Wr_ur_tang ;
	double Wr_ur_tang_z;
	
	Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
	Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
	Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
	Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
	
	//solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
	Euler3D_Polytropic_pState Wl_rotated, Wr_rotated;
	Euler3D_Polytropic_cState Flux, Flux_rotated;
	
	
	/* Apply the frame rotation and evaluate left and right
		solution states in the local rotated frame defined
		by the unit normal vector. */
	Wl_rotated.Copy(Wl);
	Wr_rotated.Copy(Wr);
	
	// Left state velocity in rotated frame
	Wl_ur_norm = dot(Wl.v, norm_dir);
	Wl_ur_tang = abs(Wl.v - Wl_ur_norm*norm_dir);
	Wl_ur_tang_vector = (Wl.v - Wl_ur_norm*norm_dir);
	if(Wl_ur_tang != ZERO){
		Wl_ur_tang_unit_vector =  Wl_ur_tang_vector/Wl_ur_tang;
	}else{
		Wl_ur_tang_unit_vector= Vector3D_ZERO;
	}
	
	Wl_rotated.rho = Wl.rho;
	Wl_rotated.v.x = Wl_ur_norm ;
	Wl_rotated.v.y = Wl_ur_tang;
	Wl_rotated.v.z = ZERO;
	Wl_rotated.p = Wl.p;
	// Right state velocity in rotated frame
	Wr_ur_norm = dot(Wr.v, norm_dir);
	Wr_ur_tang_vector = Wr.v - Wr_ur_norm*norm_dir;
	Wr_ur_tang = abs(Wr.v - Wr_ur_norm*norm_dir);
	if( Wr_ur_tang != ZERO){
		Wr_ur_tang_unit_vector =  Wr_ur_tang_vector/Wr_ur_tang ;
	}else{
		Wr_ur_tang_unit_vector= Vector3D_ZERO;  
	}
	
	Wr_rotated.rho = Wr.rho;
	Wr_rotated.v.x = Wr_ur_norm;
	Wr_rotated.v.y = dot( Wr_ur_tang_vector, Wl_ur_tang_unit_vector);
	Wr_rotated.v.z = abs( Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector);
	Wr_rotated.p = Wr.p;
	
	
	Wr_ur_tang_z = abs(Wr_ur_tang_vector-Wr_rotated.v.y* Wl_ur_tang_unit_vector);
	Wr_ur_tang_z_vector = Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector;
	if(Wr_ur_tang_z !=ZERO){
		Wr_ur_tang_z_unit_vector = Wr_ur_tang_z_vector /Wr_ur_tang_z ;
	}else{
		Wr_ur_tang_z_unit_vector = Vector3D_ZERO;
		
	}
	
	
    /* Evaluate the intermediate state solution 
		flux in the rotated frame. */
	
	Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated);
	
	/* Rotate back to the original Cartesian reference
		frame and return the solution flux. */
	
	Flux.Copy(Flux_rotated);
	
	
	
	Flux_rotated_x = Flux.rhov.x*norm_dir;
	Flux_rotated_tang_y = Flux.rhov.y* Wl_ur_tang_unit_vector ;
	Flux_rotated_tang_z = Flux.rhov.z* Wr_ur_tang_z_unit_vector;
	
	
	Flux.rhov =  Flux_rotated_x + Flux_rotated_tang_y+ Flux_rotated_tang_z;

	return (Flux);
}
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_n(const  Euler3D_Polytropic_cState &Ul,  
															   const  Euler3D_Polytropic_cState &Ur,
															   const  Vector3D &norm_dir) {
	return (FluxRoe_n(Ul.W(),Ur.W(),norm_dir));
}


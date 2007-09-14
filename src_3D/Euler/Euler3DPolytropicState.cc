/* Euler3DPolytropicState.cc:  Subroutines for 3D Euler Solution State Classes. */

/* Include 3D Euler solution state header file. */

#include "Euler3DPolytropicState.h"


/* -------------------------------------------------------------------------- *
 *				 Euler3D_Polytropic_pState subroutines						  *
 * -------------------------------------------------------------------------- */

   /*
	* Set static variables 
	* -------------------- 
	*/

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
	double Euler3D_Polytropic_pState::uo(void) const{
		return (sqrt(v.x*v.x + v.y*v.y + v.z*vz));  
	}

	// Temperature.
	double Euler3D_Polytropic_pState::T(void) {
		return (p/(rho*R));
	}
	double Euler3D_Polytropic_pState::T(void) const {
		return (p/(rho*R));
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
		return (p*gm1i + HALF*rho*sqr(uo()));
	}
	double Euler3D_Polytropic_pState::E(void) const {
		return (p*gm1i + HALF*rho*sqr(uo()));
	}

	// Specific enthalpy. 
	double Euler3D_Polytropic_pState::h(void) {
		return (g*p/(gm1*rho) + HALF*sqr(uo()));
	}
	double Euler3D_Polytropic_pState::h(void) const {
		return (g*p/(gm1*rho) + HALF*(sqr(uo())));
	}

	// Total enthalpy. 
	double Euler3D_Polytropic_pState::H(void) {
		return (g*gm1i*p + HALF*rho*sqr(uo()));
	}
	double Euler3D_Polytropic_pState::H(void) const {
		return (g*gm1i*p + HALF*rho*sqr(uo()));
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
	Vector3D Euler3D_Polytropic_pState::dv(void) {
		return (rho*v);
	}
	Vector3D Euler3D_Polytropic_pState::dv(void) const {
		return (rho*v);
		
	}
	double Euler3D_Polytropic_pState::dv(const Vector3D &n) {
		return (rho*(v*n));
		
	}
	double Euler3D_Polytropic_pState::dv(const Vector3D &n) const {
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
				*(ONE+HALF*gm1*sqr(uo())/(g*prho)));
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator +(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho+W2.rho, 
										  W1.v.x+W2.v.x, W1.v.y+W2.v.y, W1.v.z+W2.v.z,
										  W1.p+W2.p));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator -(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho-W2.rho,
										  W1.v.x-W2.v.x, W1.v.y-W2.v.y, W1.v.z-W2.v.z,
										  W1.p-W2.p));
	}
	double Euler3D_Polytropic_pState::operator *(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho*W2.rho + W1.v.x*W2.v.x + W1.v.y*W2.v.y + W1.v.z*W2.v.z + W1.p*W2.p);
	}	
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator *(const Euler3D_Polytropic_pState &W, const double &a) {
		return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z, a*W.p));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator *(const double &a, const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z,a*W.p));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator /(const Euler3D_Polytropic_pState &W, const double &a) {
		return (Euler3D_Polytropic_pState(W.rho/a, W.v.x/a, W.v.y/a, W.v.z/a,W.p/a));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator ^(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (Euler3D_Polytropic_pState(W1.rho*W2.rho,
										  W1.v.x*W2.v.x, 
										  W1.v.y*W2.v.y, 
										  W1.v.z*W2.v.z, 
										  W1.p*W2.p));
	}
	
	// Unary arithmetic operators.
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator +(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_pState(W.rho,W.v.x,W.v.y, W.v.z,W.p));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::operator -(const Euler3D_Polytropic_pState &W) {
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
		d -= W.d; 
		v.x -= W.v.x;	v.y -= W.v.y;	v.z -= W.v.z;
		p -= W.p;
		return *this;
	}
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator *=(const double &a) {
		d *= a; 
		v.x *= a;	v.y *= a	v.z *= a; 
		p *= a;
		return *this;
	}
	Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator /=(const double &a) {
		d /= a; 
		v.x /= a;	v.y /= a;	v.z /= a;
		p /= a;
		return *this;
	}

	// Relational operators.
	int Euler3D_Polytropic_pState::operator ==(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho == W2.rho && W1.v.x == W2.v.x && W1.v.y == W2.v.y
				&& W1.v.z == W2.v.z && W1.p == W2.p);
	}

	int Euler3D_Polytropic_pState::operator !=(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2) {
		return (W1.rho != W2.rho || W1.v.x != W2.v.x || W1.v.y != W2.v.y
				|| W1.v.z != W2.v.z || W1.p != W2.p);
	}

	// Input-output operators.
	ostream& Euler3D_Polytropic_pState::operator << (ostream &out_file, const Euler3D_Polytropic_pState &W) {
		out_file.setf(ios::scientific);
		out_file << " " << W.rho<< " " << W.v.x << " " << W.v.y 
			<< " " << W.v.z << " " << W.p;
		out_file.unsetf(ios::scientific);
		return (out_file);
	}
	istream& Euler3D_Polytropic_pState::operator >> (istream &in_file,  Euler3D_Polytropic_pState &W) {
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
		return (Euler3D_Polytropic_cState(rho, dv(), E()));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::U(void) const {
		return (Euler3D_Polytropic_cState(rho, dv(), E()));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::U(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho, W.rho*W.v.x, W.rho*W.v.y,
									  W.rho*W.v.z, W.E()));
	}

   /* 
	* Fluxes and Jacobians
	* --------------------
	*/

	// x-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::F(void) {
		return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::F(void) const {
		return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
	}
	Euler3D_Polytropic_cState F(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho*W.v.x, W.rho*sqr(W.v.x) + W.p, W.rho*W.v.x*W.v.y, W.rho*W.v.x*W.v.z, W.v.x*W.H()));
	}
	void Euler3D_Polytropic_pState::dFdU(DenseMatrix &dFdU) {
		dFdU(0,1) += ONE;
		dFdU(1,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.x));
		dFdU(1,1) -= v.x*(g-THREE);
		dFdU(1,2) -= v.y*gm1;
		dFdU(1,3) -= v.z*gm1;
		dFdU(1,4) += gm1;
		dFdU(2,0) -= v.x*v.y;
		dFdU(2,1) += v.y;
		dFdU(2,2) += v.x;
		dFdU(3,0) -= v.x*v.z;
		dFdU(3,1) += v.z;
		dFdU(3,3) += v.x;
		dFdU(4,0) += HALF*v.x*v.sqr()*(g-TWO) - v.x*p/rho*g*gm1i;
		dFdU(4,1) += HALF*(v.sqr() + TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
		dFdU(4,2) -= v.x*v.y*gm1;
		dFdU(4,3) -= v.x*v.z*gm1;
		dFdU(4,4) += g*v.x;
	}
	void Euler3D_Polytropic_pState::dFdU(DenseMatrix &dFdU) const {
		dFdU(0,1) += ONE;
		dFdU(1,0) += HALF*(v.sqr()*gm1 - TWO*sqr(v.x));
		dFdU(1,1) -= v.x*(g-THREE);
		dFdU(1,2) -= v.y*gm1;
		dFdU(1,3) -= v.z*gm1;
		dFdU(1,4) += gm1;
		dFdU(2,0) -= v.x*v.y;
		dFdU(2,1) += v.y;
		dFdU(2,2) += v.x;
		dFdU(3,0) -= v.x*v.z;
		dFdU(3,1) += v.z;
		dFdU(3,3) += v.x;
		dFdU(4,0) += HALF*v.x*v.sqr()*(g-TWO) - v.x*p/rho*g*gm1i;
		dFdU(4,1) += HALF*(v.sqr() + TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
		dFdU(4,2) -= v.x*v.y*gm1;
		dFdU(4,3) -= v.x*v.z*gm1;
		dFdU(4,4) += g*v.x;
	}
	void Euler3D_Polytropic_pState::dFdU(DenseMatrix &dFdU, const Euler3D_Polytropic_pState &W) {
		dFdU(0,1) += ONE;
		dFdU(1,0) += HALF*(W.v.sqr()*W.gm1 - TWO*sqr(W.v.x));
		dFdU(1,1) -= W.v.x*(W.g-THREE);
		dFdU(1,2) -= W.v.y*W.gm1;
		dFdU(1,3) -= W.v.z*W.gm1;
		dFdU(1,4) += W.gm1;
		dFdU(2,0) -= W.v.x*W.v.y;
		dFdU(2,1) += W.v.y;
		dFdU(2,2) += W.v.x;
		dFdU(3,0) -= W.v.x*W.v.z;
		dFdU(3,1) += W.v.z;
		dFdU(3,3) += W.v.x;
		dFdU(4,0) += HALF*W.v.x*W.v.sqr()*(W.g-TWO) - W.v.x*W.p/W.rho*W.g*W.gm1i;
		dFdU(4,1) += HALF*(W.v.sqr() + TWO*sqr(W.v.x)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
		dFdU(4,2) -= W.v.x*W.v.y*W.gm1;
		dFdU(4,3) -= W.v.x*W.v.z*W.gm1;
		dFdU(4,4) += W.g*W.v.x;
	}

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
		dFxdU(4,1) += HALF*(v.sqr() + TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
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
		dFxdU(4,1) += HALF*(v.sqr() + TWO*sqr(v.x)*gm1) + p/rho*g*gm1i;
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
		dFxdU(4,1) += HALF*(W.v.sqr() + TWO*sqr(W.v.x)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
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
		dFydU(3,0) -= v.x*v.z;
		dFydU(3,2) += v.z;
		dFydU(3,3) += v.y;
		dFydU(4,0) += HALF*v.y*v.sqr()*(g-TWO) - v.y*p/rho*g*gm1i;
		dFydU(4,1) -= v.x*v.y*gm1;
		dFydU(4,2) += HALF*(v.sqr() + TWO*sqr(v.y)*gm1) + p/rho*g*gm1i;
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
		dFydU(3,0) -= v.x*v.z;
		dFydU(3,2) += v.z;
		dFydU(3,3) += v.y;
		dFydU(4,0) += HALF*v.y*v.sqr()*(g-TWO) - v.y*p/rho*g*gm1i;
		dFydU(4,1) -= v.x*v.y*gm1;
		dFydU(4,2) += HALF*(v.sqr() + TWO*sqr(v.y)*gm1) + p/rho*g*gm1i;
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
		dFydU(3,0) -= W.v.x*W.v.z;
		dFydU(3,2) += W.v.z;
		dFydU(3,3) += W.v.y;
		dFydU(4,0) += HALF*W.v.y*W.v.sqr()*(W.g-TWO) - W.v.y*W.p/W.rho*W.g*W.gm1i;
		dFydU(4,1) -= W.v.x*W.v.y*W.gm1;
		dFydU(4,2) += HALF*(W.v.sqr() + TWO*sqr(W.v.y)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
		dFydU(4,3) -= W.v.y*W.v.z*W.gm1;
		dFydU(4,4) += W.g*W.v.y;
	}

	// z-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) {
		return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z + p, rho*sqr(v.z) + p, v.z*H()));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) const {
		return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z + p, rho*sqr(v.z) + p, v.z*H()));
	}
	Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W) {
		return (Euler3D_Polytropic_cState(W.rho*W.v.z, W.rho*W.v.x*W.v.z, W.rho*W.v.y*W.v.z + W.p, W.rho*sqr(W.v.z) + W.p, W.v.z*W.H()));
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
		dFydU(4,0) += HALF*v.z*v.sqr()*(g-TWO) - v.z*p/rho*g*gm1i;
		dFzdU(4,1) -= v.x*v.z*gm1;
		dFzdU(4,2) -= v.y*v.z*gm1;
		dFzdU(4,3) += HALF*(v.sqr() + TWO*sqr(v.z)*gm1) + p/rho*g*gm1i;
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
		dFydU(4,0) += HALF*v.z*v.sqr()*(g-TWO) - v.z*p/rho*g*gm1i;
		dFzdU(4,1) -= v.x*v.z*gm1;
		dFzdU(4,2) -= v.y*v.z*gm1;
		dFzdU(4,3) += HALF*(v.sqr() + TWO*sqr(v.z)*gm1) + p/rho*g*gm1i;
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
		dFydU(4,0) += HALF*W.v.z*W.v.sqr()*(W.g-TWO) - W.v.z*W.p/W.rho*W.g*W.gm1i;
		dFzdU(4,1) -= W.v.x*W.v.z*W.gm1;
		dFzdU(4,2) -= W.v.y*W.v.z*W.gm1;
		dFzdU(4,3) += HALF*(W.v.sqr() + TWO*sqr(W.v.z)*W.gm1) + W.p/W.rho*W.g*W.gm1i;
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
	void Euler3D_Polytropic_pState::dWdU(DenseMatrix &ddWdU) const {
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda(void) {
		double c = a();
		return (Euler3D_Polytropic_pState(v.x - c, v.x, v.x, v.x, v.x + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda(void) const {
		double c = a();
		return (Euler3D_Polytropic_pState(v.x - c, v.x, v.x, v.x, v.x + c));
	}
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lambda(const Euler3D_Polytropic_pState &W) {
		double c = W.a();
		return (Euler3D_Polytropic_pState(W.v.x - c, W.v.x, W.v.x, W.v.x, W.v.x + c));
	}
	double Euler3D_Polytropic_pState::lambda(int index) {
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
		};
	}
	double Euler3D_Polytropic_pState::lambda(int index) const {
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
		};
	}

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
		};
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
		};
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
		};
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
		};
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
		};
	}



   /*
	* Conserved right eigenvectors
	* ----------------------------
	*/

	// x-direction
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc(const int &index){
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, h()-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, h()+v.x*a()));
		}		
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc(const int &index) const{
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, h()-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, h()+v.x*a()));
		}		
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_x(const int &index){
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, h()-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, h()+v.x*a()));
			}		
		}
	Euler3D_Polytropic_cState Euler3D_Polytropic_pState::rc_x(const int &index) const{
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_cState(ONE, v.x-a(), v.y, v.z, h()-v.x*a()));
			case 2:
				return (Euler3D_Polytropic_cState(ONE, v.x, v.y, v.z, HALF*v.sqr())); 
			case 3:
				return (Euler3D_Polytropic_cState(ZERO, ZERO, rho, ZERO, rho*v.y));
			case 4: 
				return (Euler3D_Polytropic_cState(ZERO, ZERO, ZERO, rho, rho*v.z));
			case 5: 
				return (Euler3D_Polytropic_cState(ONE, v.x+a(), v.y, v.z, h()+v.x*a()));
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp(const int &index){
		switch(index) {
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp(const int &index) const{
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp_y(const int &index) const {
		switch(index){
			case 1: 
				return (Euler3D_Polytropic_pState(ZERO, ZERO, -rho/(TWO*a()), ZERO, HALF/a2()));
			case 2:
				return (Euler3D_Polytropic_pState(ONE, ZERO, ZERO, ZERO, -ONE/a2())); 
			case 3:
				return (Euler3D_Polytropic_pState(ZERO, ONE, ZERO, ZERO, ZERO));
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp(const int &index) {
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
	Euler3D_Polytropic_pState Euler3D_Polytropic_pState::rp(const int &index) const {
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


/* ------------------------------------------------------------------------- *
 *				 Euler3D_Polytropic_cState subroutines						 *
 * ------------------------------------------------------------------------- */

   /*
	* Set static variables 
	* -------------------- 
	*/

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
    return (R*gm1i*log(gm1*(E - HALF*rho.sqr()/rho)/pow(rho, g)));
}
double Euler3D_Polytropic_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*rho.sqr()/rho)/pow(rho, g)));
}

// Stagnation temperature. 
double Euler3D_Polytropic_cState::To(void) {
    return ((gm1*(E - HALF*rho.sqr()/rho)/(rho*R))*
			(ONE+HALF*gm1*rho.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rho.sqr()/sqr(rho)))));
}
double Euler3D_Polytropic_cState::To(void) const {
    return ((gm1*(E - HALF*rho.sqr()/rho)/(rho*R))*
			(ONE+HALF*gm1*rho.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rho.sqr()/sqr(rho)))));
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
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(void) {
		return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(void) const {
		return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rhov.x, sqr(U.rhov.x)/U.rho + U.p(), U.rhov.x*U.rhov.y/U.rho, U.rhov.x*U.rhov.z/U.rho, U.rhov.x*U.H()/U.rho));
	}
	void Euler3D_Polytropic_cState::dFdU(DenseMatrix &dFdU) {
		W().dFdU(dFdU);
	}
	void Euler3D_Polytropic_cState::dFdU(DenseMatrix &dFdU) const {
		W().dFdU(dFdU);
	}
	void Euler3D_Polytropic_cState::dFdU(DenseMatrix &dFdU, const Euler3D_Polytropic_pState &U) {
		U.W().dFdU(dFdU);
	}

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
	void Euler3D_Polytropic_cState::dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_pState &U) {
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
	void Euler3D_Polytropic_cState::dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_pState &U) {
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
	void Euler3D_Polytropic_cState::dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_pState &U) {
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
	double &operator[](int index) {
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
	const double &operator[](int index) const {
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
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator +(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho+U2.rho,U1.rhov+U2.rhov,U1.E+U2.E));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator -(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho-U2.rho,U1.rhov-U2.rhov,U1.E-U2.E));
	}

	double Euler3D_Polytropic_cState::operator *(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (U1.rho*U2.rho+U1.rhov*U2.rhov+U1.E*U2.E);
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator *(const Euler3D_Polytropic_cState &U, const double &a) {
		return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator *(const double &a, const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator /(const Euler3D_Polytropic_cState &U, const double &a) {
		return (Euler3D_Polytropic_cState(U.rho/a,U.rhov/a,U.E/a));
	}

	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator ^(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2) {
		return (Euler3D_Polytropic_cState(U1.rho*U2.rho,
										  U1.rhov.x*U2.rhov.x,
										  U1.rhov.y*U2.rhov.y,
										  U1.rhov.z*U2.rhov.z,
										  U1.E*U2.E));
	}

	// Unary arithmetic operators.
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator +(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(U.rho,U.rhov,U.E));
	}
	Euler3D_Polytropic_cState Euler3D_Polytropic_cState::operator -(const Euler3D_Polytropic_cState &U) {
		return (Euler3D_Polytropic_cState(-U.rho,-U.rhov,-U.E));
	}

	// Shortcut arithmetic operators.
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator +=(Euler3D_Polytropic_cState &U) {
		rho += U.rho;
		rhov += U.rhov;
		E += U.E;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator -=(Euler3D_Polytropic_cState &U) {
		rho -= U.rho;
		rhov -= U.rhov;
		E -= U.E;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator *=(Euler3D_Polytropic_cState &U) {
		rho *= U.rho;
		rhov *= U.rhov;
		E *= U.E;
		return *this;
	}
	Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator /=(Euler3D_Polytropic_cState &U) {
		rho /= U.rho;
		rhov /= U.rhov;
		E /= U.E;
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
	ostream& Euler3D_Polytropic_cState::operator << (ostream &out_file, const Euler3D_Polytropic_cState &U) {
		out_file.setf(ios::scientific);
		out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " "
			<< U.rhov.z << " " << U.E;
		out_file.unsetf(ios::scientific);
		return (out_file);
	}
	istream& Euler3D_Polytropic_cState::operator >> (istream &in_file, Euler3D_Polytropic_cState &U) {
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
Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_pState &Wl, const Euler3D_Polytropic_pState &Wr) {

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

Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_cState &Ul, const Euler3D_Polytropic_cState &Ur) {
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
Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_pState &Wl,
	      	          const Euler3D_Polytropic_pState &Wr) {

    double wavespeed_l, wavespeed_r;
   
    Euler3D_Polytropic_pState Wa(0), lambdas_l(0), lambdas_r(0), lambdas_a(0);
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
        Flux = Wl.F();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.F();
    } else {
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*dUrl)/
                 (wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */

    return (Flux);
}

Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
	      	          const Euler3D_Polytropic_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
}


Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_pState &Wl,
									 const Euler3D_Polytropic_pState &Wr) {
	
    double wavespeed_l, wavespeed_r;
	
    Euler3D_Polytropic_pState Wa(0), lambdas_l(0), lambdas_r(0), lambdas_a(0);
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
Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_cState &Ul,
									 const Euler3D_Polytropic_cState &Ur) {
	return (FluxHLLE_y(Ul.W(), Ur.W()));
}

Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_pState &Wl,
									 const Euler3D_Polytropic_pState &Wr) {
	
    double wavespeed_l, wavespeed_r;
	
    Euler3D_Polytropic_pState Wa(0), lambdas_l(0), lambdas_r(0), lambdas_a(0);
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

Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_cState &Ul,
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
Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_pState &Wl,
	      	          const Euler3D_Polytropic_pState &Wr,
                          const Vector3D &norm_dir) {

   double sin_beta, cos_beta, sin_alpha, cos_alpha;
   
  //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   Euler3D_Polytropic_pState Wl_rotated(), Wr_rotated();
   Euler3D_Polytropic_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

   sin_beta = norm_dir.z;
   cos_beta = sqrt(sqr(norm_dir.x)+sqr(norm_dir.y));
   if(cos_beta == ZERO){
      sin_alpha = ZERO;
      cos_alpha = ZERO;
   }else{
      sin_alpha = norm_dir.y/cos_beta;
      cos_alpha = norm_dir.x/cos_beta;
   }
     
   /* Apply the frame rotation and evaluate left and right
      solution states in the local rotated frame defined
      by the unit normal vector. */

   Wl_rotated.rho = Wl.rho;
   Wl_rotated.v.x = (Wl.v.x*cos_alpha + Wl.v.y*sin_alpha)*cos_beta+ Wl.v.z*sin_beta;
   Wl_rotated.v.y =  (Wl.v.x*cos_alpha + Wl.v.y*sin_alpha)*sin_beta+ Wl.v.z*cos_beta;
   Wl_rotated.v.z = ZERO;
   Wl_rotated.p   = Wl.p;
    
   Wr_rotated.rho = Wr.rho;
   Wr_rotated.v.x = (Wr.v.x*cos_alpha + Wr.v.y*sin_alpha)*cos_beta+ Wr.v.z*sin_beta;
   Wr_rotated.v.y =  (Wr.v.x*cos_alpha + Wr.v.y*sin_alpha)*sin_beta+ Wr.v.z*cos_beta;
   Wr_rotated.v.z = ZERO;
   Wr_rotated.p   = Wr.p;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */

    Flux.rho = Flux_rotated.d;
    Flux.rhov.x = (Flux_rotated.rhov.x*cos_beta + Flux_rotated.rhov.y*sin_beta)*cos_alpha;
    Flux.rhov.y = (Flux_rotated.rhov.x*cos_beta + Flux_rotated.rhov.y*sin_beta)*sin_alpha;
    Flux.rhov.z = (Flux_rotated.rhov.x*sin_beta + Flux_rotated.rhov.y* cos_beta);
	Flux.E = Flux_rotated.E;

    return (Flux);

}

Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
	      	          const Euler3D_Polytropic_cState &Ur,
                          const Vector3D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
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
Euler3D_Polytropic_pState Reflect(const Euler3D_Polytropic_pState &W,
	      	       const Vector3D &norm_dir) {
   
   double sin_beta, cos_beta, sin_alpha, cos_alpha;
   double dr, ur, vr, wr, pr, u, v, w;
   
   
   sin_beta = norm_dir.z;
   cos_beta = sqrt(sqr(norm_dir.x)+sqr(norm_dir.y));
   if(cos_beta == ZERO){
      sin_alpha = ZERO;
      cos_alpha = ZERO;
   }else{
      sin_alpha = norm_dir.y/cos_beta;
      cos_alpha = norm_dir.x/cos_beta;
   }
     
   /* Apply the frame rotation and calculate the primitive
      solution state variables in the local rotated frame
      defined by the unit normal vector. */

   dr = W.rho;
   ur = (W.v.x*cos_alpha + W.v.y*sin_alpha)*cos_beta+ W.v.z*sin_beta;
   vr = (W.v.x*cos_alpha + W.v.y*sin_alpha)*sin_beta+ W.v.z*cos_beta;
   wr = ZERO;
   pr = W.p;
   
   /* Reflect the normal velocity in the rotated frame. */
   ur = -ur;
   
   /* Rotate back to the original Cartesian reference frame. */

    u = (ur*cos_beta + vr*sin_beta)*cos_alpha;
    v = (ur*cos_beta + vr*sin_beta)*sin_alpha;
    w = (ur*sin_beta + vr*cos_beta);
    
    /* Return the reflected state. */

    return (Euler3D_Polytropic_pState(dr, u, v, w, pr));
       
}

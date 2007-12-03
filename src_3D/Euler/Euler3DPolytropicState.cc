/*! \file  Euler3DPolytropicState.cc
\brief Definition of member functions for 3D Euler solution state classes
       for a polytropic gas
*/

/* Include 3D Euler solution state header file. */

#include "Euler3DPolytropicState.h"

/*--------------------------------------------------------------------------------*
 *			     Euler3D_Polytropic_pState subroutines                *
 *--------------------------------------------------------------------------------*/

/*
* Set static variables
* ---------------------------- 
*/
int Euler3D_Polytropic_pState::num_vars = NUM_VAR_EULER3D;
double Euler3D_Polytropic_pState::g = GAMMA_AIR;
double Euler3D_Polytropic_pState::gm1 = GAMMA_AIR-ONE;
double Euler3D_Polytropic_pState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler3D_Polytropic_pState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
char* Euler3D_Polytropic_pState::gas_type = "AIR";

// Set gas constants.
/*!
 * Sets gas type and corresponding static variables to "AIR"
 */
void Euler3D_Polytropic_pState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
    gas_type="AIR";
}

/*!
 * Sets gas type and corresponding static variables to the specified gas
 * \param[in] string_ptr name of the gas (e.g. "AIR", "A", "CO", "CO2", 
 *                       "CH4", "H", "H2", "HE", "H2O", "N2", "O", "O2", "e")
 */
void Euler3D_Polytropic_pState::setgas(char* string_ptr) {
    gas_type=string_ptr;
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

/*
 * State functions
 * ---------------
 */
// Total velocity.
/*!
 * \f$ u0 = \sqrt{u^2 + v^2 + w^2} \f$
 */
double Euler3D_Polytropic_pState::uo(void) {
    return (sqrt(v.x*v.x + v.y*v.y + v.z*v.z));  
}

double Euler3D_Polytropic_pState::uo(void) const{
    return (sqrt(v.x*v.x + v.y*v.y + v.z*v.z));  
}

// Temperature.
/*!
 * \f$ T = \frac{p}{\rho R} \f$
 */
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
/*!
 * \f$ e = \frac{p}{(\gamma - 1) \rho} \f$
 */
double Euler3D_Polytropic_pState::e(void) {
    return (p/(gm1*rho));
}

double Euler3D_Polytropic_pState::e(void) const {
    return (p/(gm1*rho));
}

// Total energy.
/*!
 * \f$ E = \frac{p}{\gamma - 1} + \frac{1}{2} \rho \boldmath{\bar{v}}^2 \f$
 */
double Euler3D_Polytropic_pState::E(void) {
    return (p*gm1i + HALF*rho*v.sqr());
}

double Euler3D_Polytropic_pState::E(void) const {
    return (p*gm1i + HALF*rho*v.sqr());
}

// Specific enthalpy. 
/*!
 * \f$ h = \frac{p \gamma}{\rho (\gamma - 1)} + \frac{1}{2}\boldmath{\bar{v}}^2 \f$
 */
double Euler3D_Polytropic_pState::h(void) {
    return (g*gm1i*p/rho + HALF*v.sqr());
}

double Euler3D_Polytropic_pState::h(void) const {
    return (g*gm1i*p/rho + HALF*v.sqr());
}

// Total enthalpy. 
/*!
 * \f$ H = \frac{p \gamma}{\gamma - 1} + \frac{1}{2} \rho \boldmath{\bar{v}}^2 \f$
 */
double Euler3D_Polytropic_pState::H(void) {
    return (g*gm1i*p + HALF*rho*v.sqr());
}

double Euler3D_Polytropic_pState::H(void) const {
    return (g*gm1i*p + HALF*rho*v.sqr());
}

// Sound speed. 
/*!
 * \f$ a = \sqrt{\frac{\gamma p}{\rho}} \f$
 */
double Euler3D_Polytropic_pState::a(void) {
    return (sqrt(g*p/rho));
}

double Euler3D_Polytropic_pState::a(void) const {
    return (sqrt(g*p/rho));
}

// Sound speed squared. 
/*!
 * \f$ a^2 = \frac{\gamma p}{\rho} \f$
 */
double Euler3D_Polytropic_pState::a2(void) {
    return (g*p/rho);
}

double Euler3D_Polytropic_pState::a2(void) const {
    return (g*p/rho);
}

// Mach number. 
/*!
 * \f$ M = \frac{u0}{\sqrt{\frac{\gamma p}{\rho}}} \f$
 */
double Euler3D_Polytropic_pState::M(void) {
    return (uo()/sqrt(g*p/rho));
}

double Euler3D_Polytropic_pState::M(void) const {
    return (uo()/sqrt(g*p/rho));
}

// Specific entropy. 
/*!
 * \f$ s = \frac{R}{\gamma-1} \log{\frac{p}{\rho^\gamma}} \f$
 */
double Euler3D_Polytropic_pState::s(void) {
    return (R*gm1i*log(p/pow(rho, g)));
}

double Euler3D_Polytropic_pState::s(void) const {
    return (R*gm1i*log(p/pow(rho, g)));
}

// Momentum. 
/*!
 * \f$ \rho \boldmath{\bar{v}} \f$
 */
Vector3D Euler3D_Polytropic_pState::rhov(void) {
    return (rho*v);
}

Vector3D Euler3D_Polytropic_pState::rhov(void) const {
    return (rho*v);
    
}
/*!
 * \f$ \rho (\boldmath{\bar{v}}\: \boldmath{\bar{n}}) \f$
 */
double Euler3D_Polytropic_pState::rhov(const Vector3D &n) {
    return (rho*(v*n));    
}

double Euler3D_Polytropic_pState::rhov(const Vector3D &n) const {
    return (rho*(v*n));    
}

// Stagnation temperature. 
/*!
 * \f$ T0 = \frac{p}{\rho R} \left(1 + \frac{1}{2} \boldmath{\bar{v}}^2 \frac{(\gamma-1) \rho}{\gamma p}\right)   \f$
 */
double Euler3D_Polytropic_pState::To(void) {
    return ((p/(rho*R))*(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
}

double Euler3D_Polytropic_pState::To(void) const {
    return ((p/(rho*R))*(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
}

// Stagnation pressure. 
/*!
 * \f$ p0 = p (1 + \frac{1}{2} \boldmath{\bar{v}}^2 \frac{(\gamma-1) \rho}{\gamma p})^{(\frac{\gamma}{\gamma-1})}   \f$
 */
double Euler3D_Polytropic_pState::po(void) {
    return (p*pow(ONE+HALF*gm1*sqr(uo())/(g*p/rho), g*gm1i));
}

double Euler3D_Polytropic_pState::po(void) const {
    return (p*pow(ONE+HALF*gm1*sqr(uo())/(g*p/rho), g*gm1i));
}

// Stagnation sound speed. 
/*!
 * \f$ a0 = \sqrt{\frac{\gamma p}{\rho} \left(1 + \frac{1}{2} \boldmath{\bar{v}}^2 \frac{(\gamma-1) \rho}{\gamma p}\right)}  \f$
 */
double Euler3D_Polytropic_pState::ao(void) {
    return (sqrt((g*p/rho)*(ONE+HALF*gm1*sqr(uo())/(g*p/rho))));
}

double Euler3D_Polytropic_pState::ao(void) const {
    return (sqrt((g*p/rho)*(ONE+HALF*gm1*sqr(uo())/(g*p/rho))));
}

// Stagnation enthalpy.
/*!
 * \f$ h0 = \left(\frac{\gamma p}{\rho (\gamma-1)}+\frac{1}{2}\boldmath{\bar{v}}^2\right) \left(1 + \frac{1}{2} \boldmath{\bar{v}}^2 \frac{(\gamma-1) \rho}{\gamma p}\right)  \f$
 */
double Euler3D_Polytropic_pState::ho(void) {
    return ((g*p/(gm1*rho) + HALF*sqr(uo()))*
            (ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
}

double Euler3D_Polytropic_pState::ho(void) const {
    return ((g*p/(gm1*rho) + HALF*sqr(uo()))
            *(ONE+HALF*gm1*sqr(uo())/(g*p/rho)));
}

/* 
* Conserved solution state. 
* -------------------------
*/
/*! Compute the conserved solution state from the primitive solution state */
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
/*!
    \f[ Fx = 
        \left(
              \begin{array}{l}
              u \rho  \\
              \rho  u^2+p \\
              u v \rho  \\
              u w \rho  \\
              H u
              \end{array}
              \right)
        \f]
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fx(void) {
    return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fx(void) const {
    return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
}

Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_cState(W.rho*W.v.x, W.rho*sqr(W.v.x) + W.p, W.rho*W.v.x*W.v.y, W.rho*W.v.x*W.v.z, W.v.x*W.H()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::F(void) {
    return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::F(void) const {
    return (Euler3D_Polytropic_cState(rho*v.x, rho*sqr(v.x) + p, rho*v.x*v.y, rho*v.x*v.z, v.x*H()));
}

Euler3D_Polytropic_cState F(const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_cState(W.rho*W.v.x, W.rho*sqr(W.v.x) + W.p, W.rho*W.v.x*W.v.y, W.rho*W.v.x*W.v.z, W.v.x*W.H()));
}

/*! 
     \f[ \frac{dFx}{dU} = 
    \left(
      \begin{array}{lllll}
      0 & 1 & 0 & 0 & 0 \\
      \frac{1}{2} \boldmath{\bar{v}}^2 (\gamma -1)-u^2 & -u (\gamma -3) & v-v
      \gamma  & w-w \gamma  & \gamma -1 \\
      -u v & v & u & 0 & 0 \\
      -u w & w & 0 & u & 0 \\
      \frac{1}{2} u \: \boldmath{\bar{v}}^2 (\gamma -2)+\frac{p u \gamma }{\rho
          -\gamma  \rho } & \frac{1}{2} \left(\boldmath{\bar{v}}^2 - 2(\gamma-1) u^2\right)+\frac{p
              \gamma }{(\gamma -1) \rho } & -u v (\gamma -1) & -u w (\gamma -1) & u \gamma 
      \end{array}
      \right)
    \f]
*/
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
/*!
\f[ Fy = 
    \left(\begin{array}{l}
          v \rho  \\
          u v \rho  \\
          \rho  v^2+p \\
          v w \rho  \\
          H v
          \end{array}\right)
    \f]
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fy(void) {
    return (Euler3D_Polytropic_cState(rho*v.y, rho*v.x*v.y, rho*sqr(v.y) + p, rho*v.y*v.z, v.y*H()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fy(void) const {
    return (Euler3D_Polytropic_cState(rho*v.y, rho*v.x*v.y, rho*sqr(v.y) + p, rho*v.y*v.z, v.y*H()));
}

Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_cState(W.rho*W.v.y, W.rho*W.v.x*W.v.y, W.rho*sqr(W.v.y) + W.p, W.rho*W.v.y*W.v.z, W.v.y*W.H()));
}

/*!
 \f[ \frac{dFy}{dU} =
     \left(
           \begin{array}{lllll}
           0 & 0 & 1 & 0 & 0 \\
           -u v & v & u & 0 & 0 \\
           \frac{1}{2} \boldmath{\bar{v}}^2 (\gamma -1)-v^2 & u-u \gamma  & -v (\gamma -3) & w-w \gamma  & \gamma -1
           \\
           -v w & 0 & w & v & 0 \\
           \frac{1}{2} v \: \boldmath{\bar{v}}^2 (\gamma -2)+\frac{p v \gamma }{\rho -\gamma  \rho } & -u v (\gamma -1)
           & \frac{1}{2} \left(\boldmath{\bar{v}}^2 -2 (\gamma-1) v^2 \right)+\frac{p \gamma }{(\gamma -1) \rho } & -v w (\gamma -1)
           & v \gamma 
           \end{array}
           \right)
     \f]
 */
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
/*!
   \f[ Fz = 
    \left(\begin{array}{l}
          w \rho  \\
          u w \rho  \\
          v w \rho  \\
          \rho  w^2+p \\
          H w
          \end{array}\right)
    \f]
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) {
    return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z, rho*sqr(v.z) + p, v.z*H()));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::Fz(void) const {
    return (Euler3D_Polytropic_cState(rho*v.z, rho*v.x*v.z, rho*v.y*v.z, rho*sqr(v.z) + p, v.z*H()));
}

Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_cState(W.rho*W.v.z, W.rho*W.v.x*W.v.z, W.rho*W.v.y*W.v.z, W.rho*sqr(W.v.z) + W.p, W.v.z*W.H()));
}
/*!
    \f[ \frac{dFz}{dU} =
        \left(
              \begin{array}{lllll}
              0 & 0 & 0 & 1 & 0 \\
              -u w & w & 0 & u & 0 \\
              -v w & 0 & w & v & 0 \\
              \frac{1}{2} \boldmath{\bar{v}}^2 (\gamma -1)-w^2 & u-u \gamma  & v-v \gamma  & -w (\gamma -3) & \gamma -1
              \\
              \frac{1}{2} w \: \boldmath{\bar{v}}^2 (\gamma -2)+\frac{p w \gamma }{\rho -\gamma  \rho } & -u w (\gamma -1)
              & -v w (\gamma -1) & \frac{1}{2} \left(\boldmath{\bar{v}}^2 -2(\gamma-1) w^2\right)+\frac{p \gamma }{(\gamma -1) \rho }
              & w \gamma 
              \end{array}
              \right)
    \f]
 */
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
/*!
 \f[ \frac{dU}{dW} = \left(
      \begin{array}{lllll}
      1 & 0 & 0 & 0 & 0 \\
      u & \rho  & 0 & 0 & 0 \\
      v & 0 & \rho  & 0 & 0 \\
      w & 0 & 0 & \rho  & 0 \\
      \frac{1}{2} \boldmath{\bar{v}}^2 & u \rho  & v \rho  & w \rho  &
      \frac{1}{\gamma -1}
      \end{array}
      \right) \f]
 */
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
/*!
 \f[ \frac{dW}{dU} = 
     \left(
           \begin{array}{lllll}
           1 & 0 & 0 & 0 & 0 \\
           -\frac{u}{\rho } & \frac{1}{\rho } & 0 & 0 & 0 \\
           -\frac{v}{\rho } & 0 & \frac{1}{\rho } & 0 & 0 \\
           -\frac{w}{\rho } & 0 & 0 & \frac{1}{\rho } & 0 \\
           \frac{1}{2} \boldmath{\bar{v}}^2 (\gamma -1) & -u(\gamma-1)  & -v(\gamma-1)  &
           -w(\gamma-1)  & \gamma -1
           \end{array}
           \right)
     \f]
 */
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
/*!
* \f[ \lambda_x =
    \left( \begin{array}{l}
           u-a \\ u \\ u \\ u \\ u + a \end{array} \right)
    \f]
*/
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
    }
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
    }
}

// y-direction
/*!
* \f[ \lambda_y =
    \left( \begin{array}{l}
           v-a \\ v \\ v \\ v \\ v + a \end{array} \right)
    \f]
*/
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
/*!
 * \f[ \lambda_z =
     \left( \begin{array}{l}
            w-a \\ w \\ w \\ w \\ w + a \end{array} \right)
    \f]
 */
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
/*! Right conserved eigenvectors are the columns of
    \f[ Rc_x = 
        \left(
              \begin{array}{lllll}
              1 & 1 & 0 & 0 & 1 \\
              u-a & u & 0 & 0 & u+a \\
              v & v & \rho  & 0 & v \\
              w & w & 0 & \rho  & w \\
              h - u a & \frac{1}{2} \bar{v}^2 & v \rho
              & w \rho  & h + u a
              \end{array}
              \right)
    \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return right conservative eigenvector with given index
*/
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

// y-direction
/*! Right conserved eigenvectors are the columns of
\f[ Rc_y = 
    \left(
          \begin{array}{lllll}
          1 & 1 & 0 & 0 & 1 \\
          u & u & \rho  & 0 & u \\
          v-a & v & 0 & 0 & v+a \\
          w & w & 0 & \rho  & w \\
          h - v a & \frac{1}{2} \bar{v}^2 & u \rho
          & w \rho  & h + v a
          \end{array}
          \right)
    \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return right conservative eigenvector with given index
*/
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
/*! Right conserved eigenvectors are the columns of
   \f[ Rc_z = 
       \left(
             \begin{array}{lllll}
             1 & 1 & 0 & 0 & 1 \\
             u & u & \rho  & 0 & u \\
             v & v & 0 & \rho  & v \\
             w-a & w & 0 & 0 & w+a \\
             h - w a & \frac{1}{2} \bar{v}^2 & u \rho
             & v \rho  & h + w a
             \end{array}
             \right)
       \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return right conservative eigenvector with given index
*/
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
/*! Left primitive eigenvectors are the rows of
\f[ Lp_x = 
    \left(
          \begin{array}{lllll}
          0 & -\frac{\rho }{2 a} & 0 & 0 & \frac{1}{2 a^2} \\
          1 & 0 & 0 & 0 & -\frac{1}{a^2} \\
          0 & 0 & 1 & 0 & 0 \\
          0 & 0 & 0 & 1 & 0 \\
          0 & \frac{\rho }{2 a} & 0 & 0 & \frac{1}{2 a^2}
          \end{array}
          \right)
    \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return left primitive eigenvector with given index
*/
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

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp(const int &index) {
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

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::lp(const int &index) const {
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
/*! Left primitive eigenvectors are the rows of
\f[ Lp_y = 
    \left(
          \begin{array}{lllll}
          0 & 0 & -\frac{\rho }{2 a} & 0 & \frac{1}{2 a^2} \\
          1 & 0 & 0 & 0 & -\frac{1}{a^2} \\
          0 & 1 & 0 & 0 & 0 \\
          0 & 0 & 0 & 1 & 0 \\
          0 & 0 & \frac{\rho }{2 a} & 0 & \frac{1}{2 a^2}
          \end{array}
          \right)
    \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return left primitive eigenvector with given index
*/
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
/*! Left primitive eigenvectors are the rows of
\f[ Lp_z = 
    \left(
          \begin{array}{lllll}
          0 & 0 & 0 & -\frac{\rho }{2 a} & \frac{1}{2 a^2} \\
          1 & 0 & 0 & 0 & -\frac{1}{a^2} \\
          0 & 1 & 0 & 0 & 0 \\
          0 & 0 & 1 & 0 & 0 \\
          0 & 0 & 0 & \frac{\rho }{2 a} & \frac{1}{2 a^2}
          \end{array}
          \right)
    \f]
    \param [in] index The number of the requested eigenvector (1 to 5)
    \return left primitive eigenvector with given index
*/
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
/*! Right primitive eigenvectors are the columns of
\f[ Rp_x = 
    \left(
          \begin{array}{lllll}
          1 & 1 & 0 & 0 & 1 \\
          -\frac{a}{\rho } & 0 & 0 & 0 & \frac{a}{\rho } \\
          0 & 0 & 1 & 0 & 0 \\
          0 & 0 & 0 & 1 & 0 \\
          a^2 & 0 & 0 & 0 & a^2
          \end{array}
          \right)
    \f]
\param [in] index The number of the requested eigenvector (1 to 5)
\return right primitive eigenvector with given index
*/
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

// y-direction
/*! Right primitive eigenvectors are the columns of
\f[ Rp_y = 
    \left(
          \begin{array}{lllll}
          1 & 1 & 0 & 0 & 1 \\
          0 & 0 & 1 & 0 & 0 \\
          -\frac{a}{\rho } & 0 & 0 & 0 & \frac{a}{\rho } \\
          0 & 0 & 0 & 1 & 0 \\
          a^2 & 0 & 0 & 0 & a^2
          \end{array}
          \right)
    \f]
\param [in] index The number of the requested eigenvector (1 to 5)
\return right primitive eigenvector with given index
*/
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
/*! Right primitive eigenvectors are the columns of
\f[ Rp_z = 
    \left(
          \begin{array}{lllll}
          1 & 1 & 0 & 0 & 1 \\
          0 & 0 & 1 & 0 & 0 \\
          0 & 0 & 0 & 1 & 0 \\
          -\frac{a}{\rho } & 0 & 0 & 0 & \frac{a}{\rho } \\
          a^2 & 0 & 0 & 0 & a^2
          \end{array}
          \right)
    \f]
\param [in] index The number of the requested eigenvector (1 to 5)
\return right primitive eigenvector with given index
*/
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

/*
 * Flux Functions
 * --------------
 */
/*!
 * Routine: RoeAverage (Roe Averages)                   
 *                                                      
 * This function returns the Roe-averaged primitive     
 * solution state given left and right primitive        
 * solution variables.  See Roe (1981).    
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Roe averaged solution state
 */
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::RoeAverage(const Euler3D_Polytropic_pState &Wl, 
                                                                const Euler3D_Polytropic_pState &Wr) {

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

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::RoeAverage(const Euler3D_Polytropic_cState &Ul, 
                                                                const Euler3D_Polytropic_cState &Ur) {
    return (RoeAverage(Ul.W(), Ur.W()));
}

/*!
* Routine: FluxHLLE_x (Harten-Lax-van Leer flux         
*                      function, x-direction)           
*                                                       
* This function returns the intermediate state solution 
* flux for the x-direction given left and right         
* solution states by using the so-called Harten-Lax-    
* van Leer approximation for the fluxes.  See Harten,   
* Lax, van Leer (1983).       
*   \param [in] Wl The left solution state
*   \param [in] Wr The right solution state
*   \return Intermediate state solution flux in the x-direction
*/
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
    }

    /* Return solution flux. */

    return (Flux);

}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
                                                                const Euler3D_Polytropic_cState &Ur) {
    return (FluxHLLE_x(Ul.W(), Ur.W()));
}

/*!
 * Routine: FluxHLLE_y (Harten-Lax-van Leer flux         
 *                      function, y-direction)           
 *                                                       
 * This function returns the intermediate state solution 
 * flux for the y-direction given left and right         
 * solution states by using the so-called Harten-Lax-    
 * van Leer approximation for the fluxes.  See Harten,   
 * Lax, van Leer (1983).       
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Intermediate state solution flux in the y-direction
 */
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
    }

    /* Return solution flux. */

    return (Flux);
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_y(const Euler3D_Polytropic_cState &Ul,
                                                                const Euler3D_Polytropic_cState &Ur) {
    return (FluxHLLE_y(Ul.W(), Ur.W()));
}

/*!
 * Routine: FluxHLLE_z (Harten-Lax-van Leer flux         
 *                      function, z-direction)           
 *                                                       
 * This function returns the intermediate state solution 
 * flux for the z-direction given left and right         
 * solution states by using the so-called Harten-Lax-    
 * van Leer approximation for the fluxes.  See Harten,   
 * Lax, van Leer (1983).       
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Intermediate state solution flux in the z-direction
 */
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
    } 

    /* Return solution flux. */

    return (Flux);

}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_z(const Euler3D_Polytropic_cState &Ul,
                                                                const Euler3D_Polytropic_cState &Ur) {
    return (FluxHLLE_z(Ul.W(), Ur.W()));
}

/*!
* Routine: FluxHLLE_n (Harten-Lax-van Leer flux         
*                      function, n-direction)           
*                                                       
* This function returns the intermediate state solution
* flux for an arbitrary direction defined by a unit     
* normal vector in the direction of interest, given     
* left and right solution states.  The problem is       
* solved by first applying a frame rotation to rotate   
* the problem to a local frame aligned with the unit    
* normal vector and then by using the so-called         
* Harten-Lax-van Leer approximation to specify the      
* intermediate state fluxes in terms of the rotated     
* solution states.  See Harten, Lax, van Leer (1983).
* 
*   \param [in] Wl The left solution state
*   \param [in] Wr The right solution state
*   \param [in] norm_dir The direction of interest
*   \return Intermediate state solution flux in the direction norm_dir
*/
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


    /* Evaluate the intermediate state solution flux in the rotated frame. */
    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);

    /* Rotate back to the original Cartesian reference frame and return the solution flux. */
    Flux.Copy(Flux_rotated);

    Flux_rotated_x = Flux.rhov.x*norm_dir;
    Flux_rotated_tang_y = Flux.rhov.y* Wl_ur_tang_unit_vector ;
    Flux_rotated_tang_z = Flux.rhov.z* Wr_ur_tang_z_unit_vector;

    Flux.rhov =  Flux_rotated_x + Flux_rotated_tang_y + Flux_rotated_tang_z;
    return (Flux);
}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
                                                                const Euler3D_Polytropic_cState &Ur, 
                                                                const Vector3D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}

/*!
 * Routine: FluxRoe_x (Roe's flux function, x-direction)  
 *                                                      
 * This function returns the intermediate state solution 
 * flux for the x-direction given left and right         
 * solution states by using the "linearized" approximate 
 * Riemann solver of Roe for the two states.  See Roe    
 * (1981).  
 * 
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Intermediate state solution flux in the x-direction 
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_x(const Euler3D_Polytropic_pState &Wl,  
                                                               const Euler3D_Polytropic_pState &Wr) {

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
        wavespeeds = lambda_minus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1 ; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] < ZERO) {
                Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
            }
        } 

    } else {
        Flux = Wr.Fx();
        wavespeeds = lambda_plus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] > ZERO) {
                Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
            }
        } 
    } 

    /* Return solution flux. */    

    return (Flux);

}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_x(const Euler3D_Polytropic_cState &Ul,  
                                                               const Euler3D_Polytropic_cState &Ur) {
    return (FluxRoe_x(Ul.W(),Ur.W()));
}

/*!
 * Routine: FluxRoe_y (Roe's flux function, y-direction)  
 *                                                      
 * This function returns the intermediate state solution 
 * flux for the y-direction given left and right         
 * solution states by using the "linearized" approximate 
 * Riemann solver of Roe for the two states.  See Roe    
 * (1981).  
 * 
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Intermediate state solution flux in the y-direction 
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_y(const Euler3D_Polytropic_pState &Wl,  
                                                               const Euler3D_Polytropic_pState &Wr) {

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
        wavespeeds = lambda_minus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1 ; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] < ZERO)
                Flux += wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);
        } 

    } else {
        Flux = Wr.Fy();
        wavespeeds = lambda_plus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] > ZERO)
                Flux -= wavespeeds[i]*(Wa.lp_y(i)*dWrl)*Wa.rc_y(i);  
        } 
    } 

    /* Return solution flux. */    
    return (Flux);    

}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_y(const Euler3D_Polytropic_cState &Ul,  
                                                               const Euler3D_Polytropic_cState &Ur) {
    return (FluxRoe_y(Ul.W(),Ur.W()));
}

/*!
 * Routine: FluxRoe_z (Roe's flux function, z-direction)  
 *                                                      
 * This function returns the intermediate state solution 
 * flux for the z-direction given left and right         
 * solution states by using the "linearized" approximate 
 * Riemann solver of Roe for the two states.  See Roe    
 * (1981).  
 * 
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \return Intermediate state solution flux in the z-direction 
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_z(const Euler3D_Polytropic_pState &Wl,  
                                                               const Euler3D_Polytropic_pState &Wr){

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
        wavespeeds = lambda_minus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1 ; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] < ZERO)
                Flux += wavespeeds[i]*(Wa.lp_z(i)*dWrl)*Wa.rc_z(i);
        } 

    } else {
        Flux = Wr.Fz();
        wavespeeds = lambda_plus(lambdas_a, lambdas_l, lambdas_r);
        for (int i=1; i <= Wl.num_vars; i++) {
            if (wavespeeds[i] > ZERO)
                Flux -= wavespeeds[i]*(Wa.lp_z(i)*dWrl)*Wa.rc_z(i); 
        } 
    } 

    /* Return solution flux. */    
    return (Flux);    

}

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_z(const Euler3D_Polytropic_cState &Ul,  
                                                               const Euler3D_Polytropic_cState &Ur) {
    return (FluxRoe_y(Ul.W(),Ur.W()));
}

/*!
 * Routine: FluxRoe_n (Roe's flux function, n-direction) 
 *                                                       
 * This function returns the intermediate state solution 
 * flux for an arbitrary direction defined by a unit     
 * normal vector in the direction of interest, given     
 * left and right solution states.  The problem is       
 * solved by first applying a frame rotation to rotate   
 * the problem to a local frame aligned with the unit    
 * normal vector and then by using the "linearized"      
 * approximate Riemann solver of Roe to specify the flux
 * in terms of the rotated solution states.  See Roe     
 * (1981).      
 *   \param [in] Wl The left solution state
 *   \param [in] Wr The right solution state
 *   \param [in] norm_dir The direction of interest
 *   \return Intermediate state solution flux in the direction norm_dir
 */
Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_n(const Euler3D_Polytropic_pState &Wl,
                                                               const Euler3D_Polytropic_pState &Wr,
                                                               const Vector3D &norm_dir) {

    double Wl_ur_norm, Wl_ur_tang;
    double Wr_ur_norm, Wr_ur_tang ;
    double Wr_ur_tang_z;
    Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
    Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
    Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
    Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
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

Euler3D_Polytropic_cState Euler3D_Polytropic_pState::FluxRoe_n(const Euler3D_Polytropic_cState &Ul,  
                                                               const Euler3D_Polytropic_cState &Ur,
                                                               const Vector3D &norm_dir) {
    return (FluxRoe_n(Ul.W(),Ur.W(),norm_dir));
}

/******************************************************************************************
 * Euler3D_Polytropic_pState::lambda_plus -- Positive wave speeds determined using        *
 *                                           Harten entropy fix.                          *
 ******************************************************************************************/
/*!
* Routine: lambda_plus (positive wave speeds determined using Harten entropy fix)
*                                                      
*/
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::
lambda_plus(const Euler3D_Polytropic_pState &lambdas_a,
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

/******************************************************************************************
 * Euler3D_Polytropic_pState::lambda_minus -- Negative wave speeds determined using       *
 *                                            Harten entropy fix.                         *
 ******************************************************************************************/
/*!
* Routine: lambda_minus (negative wave speeds determined using Harten entropy fix) 
*                                                      
*/
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::
lambda_minus(const Euler3D_Polytropic_pState &lambdas_a,
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

/*!
 * Routine: Reflect                                     
 *                                                      
 * This function returns the reflected solution state   
 * in a given direction given the primitive solution    
 * variables and the unit normal vector in the          
 * direction of interest.
 *
 * \param [in] W The solution state we want to reflect
 * \param [in] norm_dir The unit normal vector in the direction of interest
 * \return The reflected solution state
 */
Euler3D_Polytropic_pState Euler3D_Polytropic_pState::Reflect(const Euler3D_Polytropic_pState &W, 
                                                             const Vector3D &norm_dir) {

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

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::MovingWall(const Euler3D_Polytropic_pState &Win,
                                                                const Euler3D_Polytropic_pState &Wout,
                                                                const Vector3D &norm_dir, 
                                                                const Vector3D &wall_velocity,
                                                                const Vector3D &pressure_gradient,
                                                                const int &TEMPERATURE_BC_FLAG) {

Euler3D_Polytropic_pState Temp;
Temp.Copy(Win);

if (wall_velocity == Vector3D_ZERO){
    Temp.v = -Win.v;
} else {
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
    } else {
        Temp.rho = Win.p/(Temp.Rtot()*Wout.T());
    }
}

return (Temp);

}

Euler3D_Polytropic_pState Euler3D_Polytropic_pState::NoSlip(const Euler3D_Polytropic_pState &Win,
                                                            const Euler3D_Polytropic_pState &Wout,
                                                            const Vector3D &norm_dir,
                                                            const Vector3D &pressure_gradient,
                                                            const int &TEMPERATURE_BC_FLAG) {

return (MovingWall(Win, 
                   Wout, 
                   norm_dir,
                   Vector3D_ZERO,
                   pressure_gradient,
                   TEMPERATURE_BC_FLAG));

}

/*------------------------------------------------------------------------------------*
 *				 Euler3D_Polytropic_cState subroutines	              *
 *------------------------------------------------------------------------------------*/

/*
* Set static variables 
* -------------------- 
*/
int Euler3D_Polytropic_cState::num_vars = NUM_VAR_EULER3D;
double Euler3D_Polytropic_cState::g = GAMMA_AIR;
double Euler3D_Polytropic_cState::gm1 = GAMMA_AIR-ONE;
double Euler3D_Polytropic_cState::gm1i = ONE/(GAMMA_AIR-ONE);
double Euler3D_Polytropic_cState::R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
char*  Euler3D_Polytropic_cState::gas_type = "AIR";

// set gas constants.
/*!
 * Sets gas type and corresponding static variables to "AIR"
 */
void Euler3D_Polytropic_cState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
    gas_type = "AIR";
}
/*!
 * Sets gas type and corresponding static variables to the specified gas
 * \param[in] string_ptr name of the gas (e.g. "AIR", "A", "CO", "CO2", 
 *                       "CH4", "H", "H2", "HE", "H2O", "N2", "O", "O2", "e")
 */
void Euler3D_Polytropic_cState::setgas(char *string_ptr) {
    gas_type=string_ptr;
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

/*
 * State functions
 * ---------------
 */
// Flow velocity. 
/*! \f[ \bar{v} = \frac{\bar{\rho v}}{\rho} \f] */
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
/*! \f[ p = (\gamma-1) (E - \frac{1}{2} \frac{\bar{\rho v}^2}{\rho}) \f] */
double Euler3D_Polytropic_cState::p(void) {
    return (gm1*(E - HALF*rhov.sqr()/rho));
}

double Euler3D_Polytropic_cState::p(void) const {
    return (gm1*(E - HALF*rhov.sqr()/rho));
}

// Temperature. 
/*! see Euler3D_Polytropic_pState::T(void) */
double Euler3D_Polytropic_cState::T(void) {
    return (gm1*(E - HALF*rhov.sqr()/rho)/(rho*R));
}

double Euler3D_Polytropic_cState::T(void) const {
    return (gm1*(E - HALF*rhov.sqr()/rho)/(rho*R));
}

// Specific internal energy. 
/*! see Euler3D_Polytropic_pState::e(void) */
double Euler3D_Polytropic_cState::e(void) {
    return (E/rho - HALF*rhov.sqr()/sqr(rho));
}

double Euler3D_Polytropic_cState::e(void) const {
    return (E/rho - HALF*rhov.sqr()/sqr(rho));
}

// Specific enthalpy. 
/*! see Euler3D_Polytropic_pState::h(void) */
double Euler3D_Polytropic_cState::h(void) {
    return (g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho));
}

double Euler3D_Polytropic_cState::h(void) const {
    return (g*E/rho - gm1*HALF*rhov.sqr()/sqr(rho));
}

// Total enthalpy. 
/*! see Euler3D_Polytropic_pState::H(void) */
double Euler3D_Polytropic_cState::H(void) {
    return (g*E - gm1*HALF*rhov.sqr()/rho);
}

double Euler3D_Polytropic_cState::H(void) const {
    return (g*E - gm1*HALF*rhov.sqr()/rho);
}

// Sound speed. 
/*! see Euler3D_Polytropic_pState::a(void) */
double Euler3D_Polytropic_cState::a(void) {
    return (sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))));
}

double Euler3D_Polytropic_cState::a(void) const {
    return (sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho))));
}

// Sound speed squared. 
/*! see Euler3D_Polytropic_pState::a2(void) */
double Euler3D_Polytropic_cState::a2(void) {
    return (g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)));
}

double Euler3D_Polytropic_cState::a2(void) const {
    return (g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)));
}

// Mach number.
/*! see Euler3D_Polytropic_pState::M(void) */
double Euler3D_Polytropic_cState::M(void) {
    return (abs(rhov)/(rho*sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

double Euler3D_Polytropic_cState::M(void) const {
    return (abs(rhov)/(rho*sqrt(g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

// Specific entropy. 
/*! see Euler3D_Polytropic_pState::s(void) */
double Euler3D_Polytropic_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*rhov.sqr()/rho)/pow(rho, g)));
}

double Euler3D_Polytropic_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*rhov.sqr()/rho)/pow(rho, g)));
}

// Stagnation temperature.
/*! see Euler3D_Polytropic_pState::To(void) */
double Euler3D_Polytropic_cState::To(void) {
    return ((gm1*(E - HALF*rhov.sqr()/rho)/(rho*R))*
        (ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

double Euler3D_Polytropic_cState::To(void) const {
    return ((gm1*(E - HALF*rhov.sqr()/rho)/(rho*R))*
        (ONE+HALF*gm1*rhov.sqr()/(rho*rho*g*gm1*(E/rho - HALF*rhov.sqr()/sqr(rho)))));
}

// Stagnation pressure. 
/*! see Euler3D_Polytropic_pState::po(void) */
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
/*! see Euler3D_Polytropic_pState::ao(void) */
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
/*! see Euler3D_Polytropic_pState::ho(void) */
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
/*! Compute the primitive solution state from the conservative solution state */
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
/*! see Euler3D_Polytropic_pState::Fx(void) */
Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(void) {
    return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(void) const {
    return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fx(const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(U.rhov.x, sqr(U.rhov.x)/U.rho + U.p(), U.rhov.x*U.rhov.y/U.rho, U.rhov.x*U.rhov.z/U.rho, U.rhov.x*U.H()/U.rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(void) {
    return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(void) const {
    return (Euler3D_Polytropic_cState(rhov.x, sqr(rhov.x)/rho + p(), rhov.x*rhov.y/rho, rhov.x*rhov.z/rho, rhov.x*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::F(const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(U.rhov.x, sqr(U.rhov.x)/U.rho + U.p(), U.rhov.x*U.rhov.y/U.rho, U.rhov.x*U.rhov.z/U.rho, U.rhov.x*U.H()/U.rho));
}

/*! see Euler3D_Polytropic_pState::dFxdU(DenseMatrix &dFxdU) */
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
/*! see Euler3D_Polytropic_pState::Fy(void) */
Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(void) {
    return (Euler3D_Polytropic_cState(rhov.y, rhov.x*rhov.y/rho, sqr(rhov.y)/rho + p(), rhov.y*rhov.z/rho, rhov.y*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(void) const {
    return (Euler3D_Polytropic_cState(rhov.y, rhov.x*rhov.y/rho, sqr(rhov.y)/rho + p(), rhov.y*rhov.z/rho, rhov.y*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fy(const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(U.rhov.y, U.rhov.x*U.rhov.y/U.rho, sqr(U.rhov.y)/U.rho + U.p(), U.rhov.y*U.rhov.z/U.rho, U.rhov.y*U.H()/U.rho));
}

/*! see Euler3D_Polytropic_pState::dFydU(DenseMatrix &dFydU) */
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
/*! see Euler3D_Polytropic_pState::Fz(void) */
Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(void) {
    return (Euler3D_Polytropic_cState(rhov.z, rhov.x*rhov.z/rho, rhov.y*rhov.z/rho, sqr(rhov.z)/rho + p(), rhov.z*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(void) const {
    return (Euler3D_Polytropic_cState(rhov.z, rhov.x*rhov.z/rho, rhov.y*rhov.z/rho, sqr(rhov.z)/rho + p(), rhov.z*H()/rho));
}

Euler3D_Polytropic_cState Euler3D_Polytropic_cState::Fz(const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(U.rhov.z, U.rhov.x*U.rhov.z/U.rho, U.rhov.y*U.rhov.z/U.rho, sqr(U.rhov.z)/U.rho + U.p(), U.rhov.z*U.H()/U.rho));
}

/*! see Euler3D_Polytropic_pState::dFzdU(DenseMatrix &dFzdU) */
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

// dU/dW
/*! see Euler3D_Polytropic_pState::dUdW(DenseMatrix &dUdW) */
void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW) {
    W().dUdW(dUdW);
}

void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW) const {
    W().dUdW(dUdW);
}

void Euler3D_Polytropic_cState::dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_cState &U) {
    U.W().dUdW(dUdW);
}

// dW/dU
/*! see Euler3D_Polytropic_pState::dWdU(DenseMatrix &dWdU) */
void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU) {
    W().dWdU(dWdU);
}

void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU) const {
    W().dWdU(dWdU);
}

void Euler3D_Polytropic_cState::dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_cState &U) {
    U.W().dWdU(dWdU);
}

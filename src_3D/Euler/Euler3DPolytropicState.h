/****************** Euler3DPolytropicState.h **************************
This class defines the state variables and constructors for the 
Euler3D Polytropic class.
***********************************************************************/

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED 
#define _EULER3D_POLYTROPIC_STATE_INCLUDED

class Euler3D_Polytropic_cState;
class Euler3D_Polytropic_pState;

// Required C++ libraries

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>


using namespace std;

// Required CFFC header files

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif //_TENSOR3D_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/* Define the classes. */

#define	NUM_VAR_EULER3D    5

class Euler3D_Polytropic_cState;

/**********************************************************************
 * Class: Euler3D_Polytropic_pState                                   *
 *                                                                    *
 * Member functions                                                   *
 *      d         -- Return density.                                  *
 *      v         -- Return flow velocity.                            *
 *      p         -- Return pressure.                                 *
 *      g         -- Return specific heat ratio.                      *
 *      gm1       -- Return g-1                                       *
 *      gm1i      -- Return 1/(g-1).                                  *
 *      R         -- Return gas constant.                             *
 *      setgas    -- Set gas constants.                               *
 *      T         -- Return temperature.                              *
 *      e         -- Return specific internal energy.                 * 
 *      E         -- Return total energy.                             *   
 *      h         -- Return specific enthalpy.                        *
 *      H         -- Return total enthalpy.                           *
 *      a         -- Return sound speed.                              *
 *      a2        -- Return sound speed square.                       *
 *      M         -- Return Mach number.                              *
 *      s         -- Return specific entropy.                         *
 *      dv        -- Return momentum.                                 *
 *      To        -- Return stagnation temperature.                   *
 *      po        -- Return stagnation pressure.                      *
 *      ao        -- Return stagnation sound speed.                   *
 *      ho        -- Return stagnation enthalpy.                      *
 *      U         -- Return conserved solution state.                 *
 *      F         -- Return x-direction solution flux.                *
 *      Fx        -- Return x-direction solution flux.                *
 *      Fy        -- Return y-direction solution flux.                *
 *      Fz        -- Return z-direction solution flux.                *  
 *      Fn        -- Return n-direction solution flux.                *
 *		dFdU	  -- Return x-direction jacobian					  *
 *		dFxdU	  -- Return x-direction jacobian					  *
 *		dFydU	  -- Return y-direction jacobian					  *
 *		dFzdU	  -- Return z-direction jacobian					  *
 *		dUdW	  -- Return solution variable jacobian				  *
 *		dWdU	  -- Return solution variable jacobian				  *
 *      lambda    -- Return x-direction eigenvalue.                   *
 *      lambda_x  -- Return x-direction eigenvalue.                   *
 *      lambda_y  -- Return y-direction eigenvalue.                   *
 *      lambda_z  -- Return z-direction eigenvalue.                   *
 *      rp        -- Return x-direction primitive right eigenvector.  *
 *      rp_x      -- Return x-direction primitive right eigenvector.  *
 *      rp_y      -- Return y-direction primitive right eigenvector.  *
 *      rp_z      -- Return z-direction primitive right eigenvector.  *
 *      rc        -- Return x-direction conserved right eigenvector.  *
 *      rc_x      -- Return x-direction conserved right eigenvector.  *
 *      rc_y      -- Return y-direction conserved right eigenvector.  *
 *      rc_z      -- Return z-direction conserved right eigenvector.  *
 *      lp        -- Return x-direction primitive left eigenvector.   *
 *      lp_x      -- Return x-direction primitive left eigenvector.   *
 *      lp_y      -- Return y-direction primitive left eigenvector.   *
 *      lp_z      -- Return z-direction primitive left eigenvector.   *
 *                                                                    *
 * Member operators                                                   *
 *      W		  -- a primitive solution state                       *
 *      c		  -- a scalar (double)                                *
 *                                                                    *
 *		W = W;                                                        *
 *		c = W[i];                                                     *
 *		W = W + W;                                                    *
 *		W = W - W;                                                    *
 *		c = W * W; (inner product)                                    *
 *		W = c * W;                                                    *
 *		W = W * c;                                                    *
 *		W = W / c;                                                    *
 *		W = W ^ W;                                                    *
 *		W = +W;                                                       *
 *		W = -W;                                                       *
 *		W += W;                                                       *
 *		W -= W;                                                       *
 *		W *= W;                                                       *
 *		W /= W;                                                       *
 *		W == W;                                                       *
 *		W != W;                                                       *
 *		cout << W; (output function)                                  *
 *		cin  >> W; (input function)                                   *
 *                                                                    *
 **********************************************************************/

/*********************************************************


Density:   rho  kg/m^3
Velocity:  v    m/s
Pressure:  p    Pa (N/m^2)

Molecular Mass:          M    kg/mol

Heat Capacity (const Pressure):  Cp  J/(kg*K)
Heat Capacity (const Volume):    Cv  J/(kg*K)
Specific Heat Ratio:             g
Specific Enthalpy:               h   J/kg
Specific Internal Energy:        e   J/kg
Total Enthalpy:                  H   J/kg 
Total Internal Energy:           E   J/kg

Viscosity:                       mu  kg/(m*s) N*s/m^2  
Thermal Conductivity:            k   N/(s*K)  W.(m*K)


***********************************************************
***********************************************************/
class Euler3D_Polytropic_pState{
	
  public: 
	static int num_vars;
	double		   rho;		// Density
	Vector3D	     v;		// Velocity
	double		     p;		// Pressure
	static double    g;		// Specific heat ratio.
	static double  gm1;		// g-1
	static double gm1i;		// 1/(g-1)
	static double    R;		// Gas constant.

   /*
	* Constructors
	* ------------
	*/
	
	// Creation constructor
	Euler3D_Polytropic_pState();
	
	// Copy constructor
	Euler3D_Polytropic_pState(const Euler3D_Polytropic_pState &W);
	
	// Assignment constructors
	Euler3D_Polytropic_pState(const double &d, 
							  const Vector3D &V, 
							  const double &pre);
	
	Euler3D_Polytropic_pState(const double &d, 
							  const double &vx, const double &vy, const double &vz, 
							  const double &pre);
	
	Euler3D_Polytropic_pState(const Euler3D_Polytropic_cState &U);
	
   /* 
	* Useful operators 
	* ---------------- 
	*/
	
	// Return the number of variables.
	int NumVar(void) {
		return NUM_VAR_EULER3D;
	}
	
	// Copy operator.
	void Copy(const Euler3D_Polytropic_pState &W) {
		rho = W.rho;	v = W.v;	p = W.p;
	}
	
	// Vacuum operator.
	void Vacuum(void) {
		rho = ZERO;		v = Vector3D_ZERO;		p = ZERO;
	}
	
	// Standard atmosphere operator.
	void Standard_Atmosphere(void) {
		rho = DENSITY_STDATM;	v.zero();	p = PRESSURE_STDATM;
	}
	
	// Check for unphysical state properties.
	int Unphysical_Properties(void) const {
		if (rho <= ZERO || p <= ZERO || E() <= ZERO)
			return 1;
		return 0;
	}
	
   /*
	* Set static variables 
	* -------------------- 
	*/
	
	// Set gas constants.
    void setgas(void);
    void setgas(char *string_ptr);

    // Total velocity.
    double uo(void) const;
    
    // Temperature.
    double T(void);
    double T(void) const;
	
	// gasconstant. (for compatibility reasons for tecplot output)
	double Rtot(void);
	double Rtot(void) const;
	
    // Specific internal energy.
    double e(void);
    double e(void) const;

    // Total energy.
    double E(void);
    double E(void) const;
    
    // Specific enthalpy. 
    double h(void);
    double h(void) const;

    // Total enthalpy. 
    double H(void);
    double H(void) const;

    // Sound speed. 
    double a(void);
    double a(void) const;

    // Sound speed squared. 
    double a2(void);
    double a2(void) const;
	
    // Mach number. 
    double M(void);
    double M(void) const;
    
    // Specific entropy. 
    double s(void);
    double s(void) const;

    // Momentum. 
    Vector3D rhov(void);
    Vector3D rhov(void) const;
    double rhov(const Vector3D &n);
    double rhov(const Vector3D &n) const;

    // Stagnation temperature. 
    double To(void);
    double To(void) const;

    // Stagnation pressure. 
    double po(void);
    double po(void) const;

    // Stagnation sound speed. 
    double ao(void);
    double ao(void) const;

    // Stagnation enthalpy. 
    double ho(void);
    double ho(void) const;
  
   /* 
    * Conserved solution state. 
    * -------------------------
    */
    Euler3D_Polytropic_cState U(void);
    Euler3D_Polytropic_cState U(void) const;
    Euler3D_Polytropic_cState U(const Euler3D_Polytropic_pState &W);
    
   /* 
	* Fluxes and Jacobians
	* --------------------
	*/
	
	// x-direction
	Euler3D_Polytropic_cState F(void);
    Euler3D_Polytropic_cState F(void) const;
    Euler3D_Polytropic_cState F(const Euler3D_Polytropic_pState &W);
	void dFdU(DenseMatrix &dFdU);
	void dFdU(DenseMatrix &dFdU) const;
	void dFdU(DenseMatrix &dFdU, const Euler3D_Polytropic_pState &W);
	
	Euler3D_Polytropic_cState Fx(void);
    Euler3D_Polytropic_cState Fx(void) const;
    Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_pState &W);
	void dFxdU(DenseMatrix &dFxdU);
	void dFxdU(DenseMatrix &dFxdU) const;
	void dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_pState &W);
	
	// y-direction
	Euler3D_Polytropic_cState Fy(void);
    Euler3D_Polytropic_cState Fy(void) const;
    Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_pState &W);
	void dFydU(DenseMatrix &dFydU);
	void dFydU(DenseMatrix &dFydU) const;
	void dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_pState &W);
	
	// z-direction
	Euler3D_Polytropic_cState Fz(void);
    Euler3D_Polytropic_cState Fz(void) const;
    Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W);
	void dFzdU(DenseMatrix &dFzdU);
	void dFzdU(DenseMatrix &dFzdU) const;
	void dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_pState &W);
	
   /* 
	* Solution variable Jacobians.
	* ----------------------------
	*/
	
	// dUdW
	void dUdW(DenseMatrix &dUdW);
	void dUdW(DenseMatrix &dUdW) const;
	void dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_pState &W);
	
	// dWdU
	void dWdU(DenseMatrix &dWdU);
	void dWdU(DenseMatrix &dWdU) const;
	void dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_pState &W);
		
   /* 
	* Eigenvalues
	* ----------- 
	*/
		
	// x-direction
	Euler3D_Polytropic_pState lambda(void);
    Euler3D_Polytropic_pState lambda(void) const;
    Euler3D_Polytropic_pState lambda(const Euler3D_Polytropic_pState &W);
	double lambda(int index);
	double lambda(int index) const;
	
    Euler3D_Polytropic_pState lambda_x(void);
    Euler3D_Polytropic_pState lambda_x(void) const;
    Euler3D_Polytropic_pState lambda_x(const Euler3D_Polytropic_pState &W);
	double lambda_x(int index);
	double lambda_x(int index) const;

	// y-direction
	Euler3D_Polytropic_pState lambda_y(void);
    Euler3D_Polytropic_pState lambda_y(void) const;
    Euler3D_Polytropic_pState lambda_y(const Euler3D_Polytropic_pState &W);
	double lambda_y(int index);
	double lambda_y(int index) const;
	
	// z-direction
	Euler3D_Polytropic_pState lambda_z(void);
    Euler3D_Polytropic_pState lambda_z(void) const;
    Euler3D_Polytropic_pState lambda_z(const Euler3D_Polytropic_pState &W);
	double lambda_z(int index);
	double lambda_z(int index) const;

   /*
	* Conserved right eigenvector
	* ---------------------------
	*/
	
	// x-direction
	Euler3D_Polytropic_cState rc(const int &index);
	Euler3D_Polytropic_cState rc(const int &index) const;
	
	Euler3D_Polytropic_cState rc_x(const int &index);
	Euler3D_Polytropic_cState rc_x(const int &index) const;
	
	// y-direction
	Euler3D_Polytropic_cState rc_y(const int &index);
	Euler3D_Polytropic_cState rc_y(const int &index) const;
	
	// z-direction
	Euler3D_Polytropic_cState rc_z(const int &index);
	Euler3D_Polytropic_cState rc_z(const int &index) const;
	
	
   /*
	* Primitive left eigenvector
	* --------------------------
	*/
	
	// x-direction
	Euler3D_Polytropic_pState lp(const int &index);
	Euler3D_Polytropic_pState lp(const int &index) const;
	
	Euler3D_Polytropic_pState lp_x(const int &index);
	Euler3D_Polytropic_pState lp_x(const int &index) const;
	
	// y-direction
	Euler3D_Polytropic_pState lp_y(const int &index);
	Euler3D_Polytropic_pState lp_y(const int &index) const;
	
	// z-direction
	Euler3D_Polytropic_pState lp_z(const int &index);
	Euler3D_Polytropic_pState lp_z(const int &index) const;
	
   /*
	* Primitive right eigenvector
	* ---------------------------
	*/
	
	// x-direction
	Euler3D_Polytropic_pState rp(const int &index);
	Euler3D_Polytropic_pState rp(const int &index) const;
	
	Euler3D_Polytropic_pState rp_x(const int &index);
	Euler3D_Polytropic_pState rp_x(const int &index) const;
	
	// y-direction
	Euler3D_Polytropic_pState rp_y(const int &index);
	Euler3D_Polytropic_pState rp_y(const int &index) const;
	
	// z-direction
	Euler3D_Polytropic_pState rp_z(const int &index);
	Euler3D_Polytropic_pState rp_z(const int &index) const;
		
	
   /*
	* Operators.
	* ----------
	*/
	
	// Index operator. 
	double &operator[](int index);
	const double &operator[](int index) const;
	
    // Binary arithmetic operators. 
    friend Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    friend Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    friend double operator	  				   *(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    friend Euler3D_Polytropic_pState operator *(const Euler3D_Polytropic_pState &W, const double &a);
    friend Euler3D_Polytropic_pState operator *(const double &a, const Euler3D_Polytropic_pState &W);
    friend Euler3D_Polytropic_pState operator /(const Euler3D_Polytropic_pState &W, const double &a);
	friend Euler3D_Polytropic_pState operator ^(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);

    // Unary arithmetic operators. 
    friend Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W);
    friend Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W);

    // Shortcut arithmetic operators. 
    Euler3D_Polytropic_pState &operator +=(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState &operator -=(const Euler3D_Polytropic_pState &W);
	Euler3D_Polytropic_pState &operator *=(const double &a);
	Euler3D_Polytropic_pState &operator /=(const double &a);

    // Relational operators. 
    friend int operator ==(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    friend int operator !=(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);

    // Input-output operators. 
    friend ostream& operator << (ostream &out_file, const Euler3D_Polytropic_pState &W);
    friend istream& operator >> (istream &in_file,  Euler3D_Polytropic_pState &W); 
	
	
	/* 
	 * Flux Functions
	 * --------------
	 */
	
	static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_pState &Wl,
												const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_cState &Ul,
												const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_pState &Wl,
												const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
												const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_pState &Wl,
												const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_cState &Ul,
												const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_pState &Wl,
												const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_cState &Ul,
												const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_pState &Wl,
												const Euler3D_Polytropic_pState &Wr,
												const Vector3D &norm_dir);
	static Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
												const Euler3D_Polytropic_cState &Ur,
												const Vector3D &norm_dir);
	
	static Euler3D_Polytropic_cState FluxRoe_x(const Euler3D_Polytropic_pState &Wl, 
											   const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxRoe_x(const Euler3D_Polytropic_cState &Ul, 
											   const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxRoe_y(const Euler3D_Polytropic_pState &Wl, 
											   const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxRoe_y(const Euler3D_Polytropic_cState &Ul, 
											   const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxRoe_z(const Euler3D_Polytropic_pState &Wl, 
											   const Euler3D_Polytropic_pState &Wr);
	static Euler3D_Polytropic_cState FluxRoe_z(const Euler3D_Polytropic_cState &Ul, 
											   const Euler3D_Polytropic_cState &Ur);
	static Euler3D_Polytropic_cState FluxRoe_n(const Euler3D_Polytropic_pState &Wl,
											   const Euler3D_Polytropic_pState &Wr,
											   const Vector3D &norm_dir);
	static Euler3D_Polytropic_cState FluxRoe_n(const Euler3D_Polytropic_cState &Ul,
											   const Euler3D_Polytropic_cState &Ur,
											   const Vector3D &norm_dir);
	
	friend Euler3D_Polytropic_pState HartenFixNeg(const Euler3D_Polytropic_pState  &lambda_a,
												  const Euler3D_Polytropic_pState  &lambda_l,
												  const Euler3D_Polytropic_pState  &lambda_r);
	friend Euler3D_Polytropic_pState HartenFixPos(const Euler3D_Polytropic_pState  &lambda_a,
												  const Euler3D_Polytropic_pState  &lambda_l,
												   const Euler3D_Polytropic_pState  &lambda_r);
	
	
	/*
	 * Boundary Conditions
	 * -------------------
	 */
	static Euler3D_Polytropic_pState Reflect(const Euler3D_Polytropic_pState &W, const Vector3D &norm_dir);
	static Euler3D_Polytropic_pState Moving_Wall(
										 const Euler3D_Polytropic_pState &Win,
										 const Euler3D_Polytropic_pState &Wout,
										 const Vector3D &norm_dir, 
										 const Vector3D &wall_velocity,
										 const Vector3D &pressure_gradient,
										 const int &TEMPERATURE_BC_FLAG);
	static Euler3D_Polytropic_pState No_Slip(
									 const Euler3D_Polytropic_pState &Win,
									 const Euler3D_Polytropic_pState &Wout,
									 const Vector3D &norm_dir,
									 const Vector3D &pressure_gradient,
									 const int &TEMPERATURE_BC_FLAG);
    
};

/********************************************************
 * Class: Euler3D_Polytropic_cState                     *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      E       -- Return total energy.                 *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1.                          *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          *
 *      W       -- Return primitive solution state.     *
 *                                                      *
 * Member operators                                     *
 *      U -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * U = U;                                               *
 * c = U[i];                                            *
 * U = U + U;                                           *
 * U = U - U;                                           *
 * c = U * U; (inner product)                           *
 * U = c * U;                                           *
 * U = U * c;                                           *
 * U = U / c;                                           *
 * U = +U;                                              *
 * U = -U;                                              *
 * U += U;                                              *
 * U -= U;                                              *
 * U == U;                                              *
 * U != U;                                              *
 * cout << U; (output function)                         *
 * cin  >> U; (input function)                          *
 *                                                      *
 ********************************************************/
class Euler3D_Polytropic_cState{
  
  public:
    double        rho;   // Density.
    Vector3D     rhov;   // Momentum.
    double          E;   // Total Energy.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
    static double   R;   // Gas constant.
	static int num_vars; // Number of variables.
	
	/*
	 * Constructors
	 * ------------
	 */
	
	// Creation constructor
	Euler3D_Polytropic_cState();
	
	// Copy constructor
	Euler3D_Polytropic_cState(const Euler3D_Polytropic_cState &U);
	
	// Assignment constructors
	Euler3D_Polytropic_cState(const double &d, 
							  const Vector3D &dv, 
							  const double &Etotal);
	
	Euler3D_Polytropic_cState(const double &d, 
							  const double &dvx, const double &dvy, const double &dvz, 
							  const double &Etotal);
	
	Euler3D_Polytropic_cState(const Euler3D_Polytropic_pState &W);
	
   /* 
	* Useful operators 
	* ---------------- 
	*/
	
	// Return the number of variables.
	int NumVar(void) {
		return NUM_VAR_EULER3D;
	}
	
	// Copy operator.
	void Copy(const Euler3D_Polytropic_cState &U) {
		rho = U.rho;	rhov = U.rhov;		E = U.E;
	}
	
	// Vacuum operator.
	void Vacuum(void) {
		rho = ZERO;		rhov = Vector3D_ZERO;		E = ZERO;
	}
	
	// Standard atmosphere operator.
	void Standard_Atmosphere(void) {
		rho = DENSITY_STDATM;	rhov.zero();	E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
	}
	
	// Check for unphysical state properties.
	int Unphysical_Properties(void) const {
		if (rho <= ZERO || E <= ZERO || e() <= ZERO)
			return 1;
		return 0;
	}
	
        
   /*
	* Set static variables 
	* -------------------- 
	*/
	
    // Set gas constants. 
    void setgas(void);
    void setgas(char *string_ptr);

    // Flow velocity. 
    Vector3D v(void);
    Vector3D v(void) const;
    double v(const Vector3D &n);
    double v(const Vector3D &n) const;
    
    // Pressure. 
    double p(void);
    double p(void) const;

    // Temperature. 
    double T(void);
    double T(void) const;

    // Specific internal energy. 
    double e(void);
    double e(void) const;

    // Specific enthalpy. 
    double h(void);
    double h(void) const;

    // Total enthalpy. 
    double H(void);
    double H(void) const;

    // Sound speed. 
    double a(void);
    double a(void) const;

    // Sound speed squared. 
    double a2(void);
    double a2(void) const;

    // Mach number. 
    double M(void);
    double M(void) const;

    // Specific entropy. 
    double s(void);
    double s(void) const;

    // Stagnation temperature. 
    double To(void);
    double To(void) const;

    // Stagnation pressure. 
    double po(void);
    double po(void) const;

    // Stagnation sound speed. 
    double ao(void);
    double ao(void) const;

    // Stagnation enthalpy. 
    double ho(void);
    double ho(void) const;

   /*
	* Primitive solution state. 
	* -------------------------
	*/
	
    Euler3D_Polytropic_pState W(void);
    Euler3D_Polytropic_pState W(void) const;
    Euler3D_Polytropic_pState W(const Euler3D_Polytropic_cState &U);
    
   /* 
	* Fluxes and Jacobians
	* --------------------
	*/
	
	// x-direction
	Euler3D_Polytropic_cState Fx(void);
    Euler3D_Polytropic_cState Fx(void) const;
    Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_cState &U);
	void dFxdU(DenseMatrix &dFxdU);
	void dFxdU(DenseMatrix &dFxdU) const;
	void dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_cState &U);
	
	// y-direction
	Euler3D_Polytropic_cState Fy(void);
    Euler3D_Polytropic_cState Fy(void) const;
    Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_cState &U);
	void dFydU(DenseMatrix &dFydU);
	void dFydU(DenseMatrix &dFydU) const;
	void dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_cState &U);
	
	// z-direction
	Euler3D_Polytropic_cState Fz(void);
    Euler3D_Polytropic_cState Fz(void) const;
    Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_cState &U);
	void dFzdU(DenseMatrix &dFzdU);
	void dFzdU(DenseMatrix &dFzdU) const;
	void dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_cState &U);
	
   /* 
	* Solution variable Jacobians.
	* ----------------------------
	*/
	
	// dUdW
	void dUdW(DenseMatrix &dUdW);
	void dUdW(DenseMatrix &dUdW) const;
	void dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_cState &U);
	
	// dWdU
	void dWdU(DenseMatrix &dWdU);
	void dWdU(DenseMatrix &dWdU) const;
	void dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_cState &U);
	
	
   /*
	* Operators.
	* ----------
	*/
	
	// Index operator. 
    double& operator[](int index);
    const double& operator[](int index) const;

    // Binary arithmetic operators.
    friend Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    friend Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    friend double operator					   *(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    friend Euler3D_Polytropic_cState operator *(const Euler3D_Polytropic_cState &U, const double &a);
    friend Euler3D_Polytropic_cState operator *(const double &a, const Euler3D_Polytropic_cState &U);
    friend Euler3D_Polytropic_cState operator /(const Euler3D_Polytropic_cState &U, const double &a);
	friend Euler3D_Polytropic_cState operator ^(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);

    // Unary arithmetic operators.
    friend Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U);
    friend Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U);

    // Shortcut arithmetic operators.
    Euler3D_Polytropic_cState& operator +=(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState& operator -=(const Euler3D_Polytropic_cState &U);
	Euler3D_Polytropic_cState& operator *=(const double &a);
	Euler3D_Polytropic_cState& operator /=(const double &a);
    
    // Relational operators.
    friend int operator ==(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    friend int operator !=(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    
    // Input-output operators.
    friend ostream& operator << (ostream &out_file, const Euler3D_Polytropic_cState &U);
    friend istream& operator >> (istream &in_file,  Euler3D_Polytropic_cState &U);
       
};




/* -------------------------------------------------------------------------- *
 *                     Useful 3D Euler state constants.                       *
 * -------------------------------------------------------------------------- */
const Euler3D_Polytropic_pState Euler3D_W_STDATM(DENSITY_STDATM,
				      Vector3D_ZERO, PRESSURE_STDATM);
const Euler3D_Polytropic_pState Euler3D_W_VACUUM(ZERO, Vector3D_ZERO, ZERO);
const Euler3D_Polytropic_cState Euler3D_U_STDATM(Euler3D_W_STDATM);
const Euler3D_Polytropic_cState Euler3D_U_VACUUM(Euler3D_W_VACUUM);

///* -------------------------------------------------------------------------- *
// *                Flux functions in Euler3DPolytropicState.cc                 *
// * -------------------------------------------------------------------------- */
//
//extern Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_pState &Wl,
//	      	                 const Euler3D_Polytropic_pState &Wr);
//extern Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_pState &Wl,
//	      	                 const Euler3D_Polytropic_pState &Wr);
//extern Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
//	      	                 const Euler3D_Polytropic_cState &Ur);
//extern Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_pState &Wl,
//	      	                 const Euler3D_Polytropic_pState &Wr,
//                                 const Vector3D &norm_dir);
//extern Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
//	      	                 const Euler3D_Polytropic_cState &Ur,
//                                 const Vector3D &norm_dir);
//extern Euler3D_Polytropic_pState Reflect(const Euler3D_Polytropic_pState &W,
//	      	              const Vector3D &norm_dir);

#endif /* _EULER3D_POLYTROPIC_STATE_INCLUDED  */

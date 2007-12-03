/*! \file Euler3DPolytropicState.h
 * 	\brief	Header file defining the Euler solution state classes 
 *              associated with solution of compressible inviscid flows 
 *              governed by the polytropic equation of state.
 */

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED 
#define _EULER3D_POLYTROPIC_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include required CFFC header files. */

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

class Euler3D_Polytropic_pState;
class Euler3D_Polytropic_cState;

/*! Number of variables in Euler3D_Polytropic state classes */
#define	NUM_VAR_EULER3D    5

/*! 
 * Class: Euler3D_Polytropic_pState
 *
 * \brief Primitive state solution class for 3D Euler equations
 *        governing polytropic gaseous flows.
 *
 *	This class contains the primitive variables and member functions.
 *	It also contains the approximate riemann solvers.
 *
 * Member functions
 *  - rho           -- Return density.
 *  - v             -- Return flow velocity.
 *  - p             -- Return pressure.
 *  - g             -- Return specific heat ratio.
 *  - gm1           -- Return g-1
 *  - gm1i          -- Return 1/(g-1).
 *  - R             -- Return gas constant.
 *  - setgas        -- Set gas constants.
 *  - T             -- Return temperature.
 *  - e             -- Return specific internal energy.
 *  - E             -- Return total energy.
 *  - h             -- Return specific enthalpy.
 *  - H             -- Return total enthalpy.
 *  - a             -- Return sound speed.
 *  - a2            -- Return sound speed square.
 *  - M             -- Return Mach number.
 *  - s             -- Return specific entropy.
 *  - rhov          -- Return momentum.
 *  - To            -- Return stagnation temperature.
 *  - po            -- Return stagnation pressure.
 *  - ao            -- Return stagnation sound speed.
 *  - ho            -- Return stagnation enthalpy.
 *  - U             -- Return conserved solution state.
 *  - Fx            -- Return x-direction solution flux.
 *  - Fy            -- Return y-direction solution flux.
 *  - Fz            -- Return z-direction solution flux.
 *  - dFxdU         -- Return x-direction jacobian.
 *  - dFydU         -- Return y-direction jacobian.
 *  - dFzdU         -- Return z-direction jacobian.
 *  - dUdW          -- Return solution variable jacobian.
 *  - dWdU          -- Return solution variable jacobian.
 *  - lambda_x      -- Return x-direction eigenvalue.
 *  - lambda_y      -- Return y-direction eigenvalue.
 *  - lambda_z      -- Return z-direction eigenvalue.
 *  - rp_x          -- Return x-direction primitive right eigenvector.
 *  - rp_y          -- Return y-direction primitive right eigenvector.
 *  - rp_z          -- Return z-direction primitive right eigenvector.
 *  - rc_x          -- Return x-direction conserved right eigenvector.
 *  - rc_y          -- Return y-direction conserved right eigenvector.
 *  - rc_z          -- Return z-direction conserved right eigenvector.
 *  - lp_x          -- Return x-direction primitive left eigenvector.
 *  - lp_y          -- Return y-direction primitive left eigenvector.
 *  - lp_z          -- Return z-direction primitive left eigenvector.
 *  - RoeAverage    -- Return Roe average.
 *  - HartenFixPos  -- Return positive entropy fix of Harten.
 *  - HartenFixNeg  -- Return negative entropy fix of Harten.
 *  - FluxRoe_x     -- Return Roe's solution of Riemann problem in x-direction.
 *  - FluxRoe_y     -- Return Roe's solution of Riemann problem in y-direction.
 *  - FluxRoe_z     -- Return Roe's solution of Riemann problem in z-direction.
 *  - FluxRoe_n     -- Return Roe's solution of Riemann problem in n-direction.
 *  - FluxHLLE_x    -- Return HLLE solution of Riemann problem in x-direction.
 *  - FluxHLLE_y    -- Return HLLE solution of Riemann problem in y-direction.
 *  - FluxHLLE_z    -- Return HLLE solution of Riemann problem in z-direction.
 *  - FluxHLLE_n    -- Return HLLE solution of Riemann problem in n-direction.
 *
 * Member operators \n
 *  W         -- a primitive solution state \n
 *  c         -- a scalar (double)
 *
 *  - W = W;
 *  - c = W[i];
 *  - W = W + W;
 *  - W = W - W;
 *  - c = W * W; (inner product)
 *  - W = c * W;
 *  - W = W * c;
 *  - W = W / c;
 *  - W = W ^ W;
 *  - W = -W;
 *  - W += W;
 *  - W -= W;
 *  - W *= W;
 *  - W /= W;
 *  - W == W;
 *  - W != W;
 *  - cout << W; (output function)
 *  - cin  >> W; (input function)
 *  \nosubgrouping
 */
class Euler3D_Polytropic_pState{
public: 
	static int   num_vars;      //!< Number of variables 
	double            rho;      //!< Density                   kg/m^3
	Vector3D            v;      //!< Velocity                  m/s
	double              p;      //!< Pressure                  Pa 
	static double       g;      //!< Specific heat ratio g 
	static double     gm1;      //!< g-1 
	static double    gm1i;      //!< 1/(g-1) 
	static double       R;      //!< Gas constant R  
	static char* gas_type;      //!< Gas type

/** @name Constructors and desctructors */
/*        ----------------------------- */
//@{
    //! Creation constructor 
    Euler3D_Polytropic_pState(void): 
      rho(DENSITY_STDATM), p(PRESSURE_STDATM) {
       v.zero(); p = PRESSURE_STDATM;
    }

    //! Copy constructor 
    Euler3D_Polytropic_pState(const Euler3D_Polytropic_pState &W):
      rho(W.rho), v(W.v), p(W.p) {
    }

    //! Assignment constructor 
    Euler3D_Polytropic_pState(const double &d, 
                              const Vector3D &V, 
                              const double &pre):
      rho(d), v(V), p(pre) {
    }

    //! Assignment constructor 
    Euler3D_Polytropic_pState(const double &d, 
                              const double &vx, 
                              const double &vy, 
                              const double &vz, 
                              const double &pre):
      rho(d), p(pre) {
      v.x=vx; v.y=vy; v.z = vz;
    }

    //! Assignment constructor 
    Euler3D_Polytropic_pState(const Euler3D_Polytropic_cState &U1);
//@}
        
/** @name Useful operators */
/*        ---------------- */
//@{
    //! Return the number of variables.
    int NumVar(void) {
        return NUM_VAR_EULER3D;
    }

    //! Copy operator.
    /*! Set this primitive state to a given primitive state */
    void Copy(const Euler3D_Polytropic_pState &W) {
        rho = W.rho;   v = W.v;   p = W.p;
    }

    //! Vacuum operator.
    /*! Set this primitive state to vacuum state */
    void Vacuum(void) {
        rho = ZERO;   v = Vector3D_ZERO;   p = ZERO;
    }

    //! Standard atmosphere operator.
    /*! Set this primitive state to the standard atmosphere */
    void Standard_Atmosphere(void) {
        rho = DENSITY_STDATM;   v.zero();   p = PRESSURE_STDATM;
    }

    //! Check for unphysical state properties.
    int Unphysical_Properties(void) const {
        if (rho <= ZERO || p <= ZERO || E() <= ZERO)
            return 1;
        return 0;
    }
//@}
        
/** @name Set static variables */
/*        -------------------- */
//@{
    void setgas(void);                //!< Set gas to air.
    void setgas(char* string_ptr);    //!< Set gas.
//@}
        
/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
    double uo(void);                    //!< Total velocity.
    double T(void);                     //!< Temperature.
    double Rtot(void);                  //!< gasconstant. (for compatibility)
    double e(void);                     //!< Specific internal energy.
    double E(void);                     //!< Total energy.
    double h(void);                     //!< Specific enthalpy.
    double H(void);                     //!< Total enthalpy.
    double a(void);                     //!< Sound speed.
    double a2(void);                    //!< Sound speed squared.
    double M(void);                     //!< Mach number.
    double s(void);                     //!< Specific entropy. 
    Vector3D rhov(void);                //!< Momentum.
    double rhov(const Vector3D &n);     //!< Momentum in given direction.
    double To(void);                    //!< Stagnation temperature.
    double po(void);                    //!< Stagnation pressure.
    double ao(void);                    //!< Stagnation sound speed.
    double ho(void);                    //!< Stagnation enthalpy. 
//@}
        
    // Same functions with const 
    double uo(void) const;
    double T(void) const;
    double Rtot(void) const;
    double e(void) const;
    double E(void) const;
    double h(void) const;
    double H(void) const;
    double a(void) const;
    double a2(void) const;
    double M(void) const;
    double s(void) const;
    Vector3D rhov(void) const;
    double rhov(const Vector3D &n) const;
    double To(void) const;
    double po(void) const;
    double ao(void) const;
    double ho(void) const;    
        
/** @name Conserved solution state. */ 
/*        ------------------------- */
//@{
    //! Convert primitive solution state to conservative solution state
    Euler3D_Polytropic_cState U(void);
    static Euler3D_Polytropic_cState U(const Euler3D_Polytropic_pState &W);
//@}

    // same functions with const
    Euler3D_Polytropic_cState U(void) const;

        
/** @name Fluxes */
/*        ------ */
//@{
    Euler3D_Polytropic_cState F(void);         //!< x-direction inviscid solution flux
    Euler3D_Polytropic_cState Fx(void);         //!< x-direction inviscid solution flux
    Euler3D_Polytropic_cState Fy(void);         //!< y-direction inviscid solution flux
    Euler3D_Polytropic_cState Fz(void);         //!< z-direction inviscid solution flux
//@}

    //same functions with const
    Euler3D_Polytropic_cState F(void) const;
    static Euler3D_Polytropic_cState F(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_cState Fx(void) const;
    static Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_cState Fy(void) const;
    static Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_cState Fz(void) const;
    static Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W);
                
/** @name Flux Jacobians */
/*        -------------- */
//@{
    void dFxdU(DenseMatrix &dFxdU);         //!< x-direction flux Jacobian
    void dFydU(DenseMatrix &dFydU);         //!< y-direction flux Jacobian
    void dFzdU(DenseMatrix &dFzdU);         //!< z-direction flux Jacobian
//@}

    //same functions with const
    void dFxdU(DenseMatrix &dFxdU) const;
    static void dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_pState &W);
    void dFydU(DenseMatrix &dFydU) const;
    static void dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_pState &W);
    void dFzdU(DenseMatrix &dFzdU) const;
    static void dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_pState &W);
        
/** @name Solution variable Jacobians. */
/*        ---------------------------- */
//@{
    void dUdW(DenseMatrix &dUdW);         //!< dU/dW
    void dWdU(DenseMatrix &dWdU);         //!< dW/dU
//@}

    // same functions with const
    void dUdW(DenseMatrix &dUdW) const;
    static void dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_pState &W);
    void dWdU(DenseMatrix &dWdU) const;
    static void dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_pState &W);
        
/** @name Eigenvalues */
/*        ----------- */
//@{	
    Euler3D_Polytropic_pState lambda(void);     //!< x-direction eigenvalues
    double lambda(int index);                   //!< x-direction eigenvalues

    Euler3D_Polytropic_pState lambda_x(void);   //!< x-direction eigenvalues
    double lambda_x(int index);                 //!< x-direction eigenvalues

    Euler3D_Polytropic_pState lambda_y(void);   //!< y-direction eigenvalues
    double lambda_y(int index);                 //!< y-direction eigenvalues

    Euler3D_Polytropic_pState lambda_z(void);   //!< z-direction eigenvalues
    double lambda_z(int index);                 //!< z-direction eigenvalues
//@}

    //same functions with const
    Euler3D_Polytropic_pState lambda(void) const;
    static Euler3D_Polytropic_pState lambda(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState lambda_x(void) const;
    static Euler3D_Polytropic_pState lambda_x(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState lambda_y(void) const;
    static Euler3D_Polytropic_pState lambda_y(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState lambda_z(void) const;
    static Euler3D_Polytropic_pState lambda_z(const Euler3D_Polytropic_pState &W);
    double lambda(int index) const;
    double lambda_x(int index) const;
    double lambda_y(int index) const;
    double lambda_z(int index) const;

/** @name Conserved right eigenvector */
/*        --------------------------- */
//@{
    Euler3D_Polytropic_cState rc(const int &index);     //!< x-direction conserved right eigenvector
    Euler3D_Polytropic_cState rc_x(const int &index);   //!< x-direction conserved right eigenvector
    Euler3D_Polytropic_cState rc_y(const int &index);   //!< y-direction conserved right eigenvector
    Euler3D_Polytropic_cState rc_z(const int &index);   //!< z-direction conserved right eigenvector
//@}

    //same functions with const
    Euler3D_Polytropic_cState rc(const int &index) const;
    Euler3D_Polytropic_cState rc_x(const int &index) const;
    Euler3D_Polytropic_cState rc_y(const int &index) const;
    Euler3D_Polytropic_cState rc_z(const int &index) const;

/** @name  Primitive left eigenvector */
/*         -------------------------- */
//@{
    Euler3D_Polytropic_pState lp(const int &index);     //!< x-direction primitive left eigenvector
    Euler3D_Polytropic_pState lp_x(const int &index);   //!< x-direction primitive left eigenvector
    Euler3D_Polytropic_pState lp_y(const int &index);   //!< y-direction primitive left eigenvector
    Euler3D_Polytropic_pState lp_z(const int &index);   //!< z-direction primitive left eigenvector
//@}

    //same functions with const
    Euler3D_Polytropic_pState lp(const int &index) const;
    Euler3D_Polytropic_pState lp_x(const int &index) const;
    Euler3D_Polytropic_pState lp_y(const int &index) const;
    Euler3D_Polytropic_pState lp_z(const int &index) const;

/** @name Primitive right eigenvector */
/*        --------------------------- */
//@{
    Euler3D_Polytropic_pState rp(const int &index);     //!< x-direction primitive right eigenvector
    Euler3D_Polytropic_pState rp_x(const int &index);   //!< x-direction primitive right eigenvector
    Euler3D_Polytropic_pState rp_y(const int &index);   //!< y-direction primitive right eigenvector
    Euler3D_Polytropic_pState rp_z(const int &index);   //!< z-direction primitive right eigenvector
//@}

    //same functions with const
    Euler3D_Polytropic_pState rp(const int &index) const;
    Euler3D_Polytropic_pState rp_x(const int &index) const;
    Euler3D_Polytropic_pState rp_y(const int &index) const;
    Euler3D_Polytropic_pState rp_z(const int &index) const;

        
/** @name Operators. */
/*        ---------- */
//@{
// Index operator. 
    //!  c = W[i]
    double &operator[](int index);		
    const double &operator[](int index) const;
        
// Binary arithmetic operators. 
    //! W = W + W
    friend Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    //! W = W - W
    friend Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    //! c = W * W  (inner product)
    friend double operator                    *(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    //! W = W * c
    friend Euler3D_Polytropic_pState operator *(const Euler3D_Polytropic_pState &W, const double &a);
    //! W = c * W
    friend Euler3D_Polytropic_pState operator *(const double &a, const Euler3D_Polytropic_pState &W);
    //! W = W / c
    friend Euler3D_Polytropic_pState operator /(const Euler3D_Polytropic_pState &W, const double &a);
    //! W = W ^ W
    friend Euler3D_Polytropic_pState operator ^(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);

// Unary arithmetic operators. 
    //! W = -W 
    friend Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W);

// Shortcut arithmetic operators. 
    //! W += W
    Euler3D_Polytropic_pState &operator +=(const Euler3D_Polytropic_pState &W);
    //! W -= W
    Euler3D_Polytropic_pState &operator -=(const Euler3D_Polytropic_pState &W);
    //! W *= W
    Euler3D_Polytropic_pState &operator *=(const double &a);
    //! W /= W
    Euler3D_Polytropic_pState &operator /=(const double &a);

// Relational operators. 
    //! W == W
    friend int operator ==(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);
    //! W != W
    friend int operator !=(const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2);

// Input-output operators. 
    //! ostream << W
    friend ostream& operator << (ostream &out_file, const Euler3D_Polytropic_pState &W);
    //! istream >> W
    friend istream& operator >> (istream &in_file,  Euler3D_Polytropic_pState &W); 
//@}
        
/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
//@{
    //! HLLE flux function in x-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr);
    //! HLLE flux function in x-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxHLLE_x(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur);
    //! HLLE flux function in y-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr);
    //! HLLE flux function in y-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxHLLE_y(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur);
    //! HLLE flux function in z-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr);
    //! HLLE flux function in z-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxHLLE_z(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur);
    //! HLLE flux function in n-direction given 2 primitive states and a direction
    static Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr,
                                                const Vector3D &norm_dir);
    //! HLLE flux function in n-direction given 2 conservative states and a direction
    static Euler3D_Polytropic_cState FluxHLLE_n(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur,
                                                const Vector3D &norm_dir);

    //! Roe flux function in x-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxRoe_x(const Euler3D_Polytropic_pState &Wl, 
                                               const Euler3D_Polytropic_pState &Wr);
    //! Roe flux function in x-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxRoe_x(const Euler3D_Polytropic_cState &Ul, 
                                               const Euler3D_Polytropic_cState &Ur);
    //! Roe flux function in y-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxRoe_y(const Euler3D_Polytropic_pState &Wl, 
                                               const Euler3D_Polytropic_pState &Wr);
    //! Roe flux function in y-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxRoe_y(const Euler3D_Polytropic_cState &Ul, 
                                               const Euler3D_Polytropic_cState &Ur);
    //! Roe flux function in y-direction given 2 primitive states
    static Euler3D_Polytropic_cState FluxRoe_z(const Euler3D_Polytropic_pState &Wl, 
                                               const Euler3D_Polytropic_pState &Wr);
    //! Roe flux function in y-direction given 2 conservative states
    static Euler3D_Polytropic_cState FluxRoe_z(const Euler3D_Polytropic_cState &Ul, 
                                               const Euler3D_Polytropic_cState &Ur);
    //! Roe flux function in n-direction given 2 primitive states and a direction
    static Euler3D_Polytropic_cState FluxRoe_n(const Euler3D_Polytropic_pState &Wl,
                                               const Euler3D_Polytropic_pState &Wr,
                                               const Vector3D &norm_dir);
    //! Roe flux function in n-direction given 2 conservative states and a direction
    static Euler3D_Polytropic_cState FluxRoe_n(const Euler3D_Polytropic_cState &Ul,
                                               const Euler3D_Polytropic_cState &Ur,
                                               const Vector3D &norm_dir);
    //@}
    //! HLLE wavespeeds in n-direction given 2 primitive states and a direction
    static Vector2D HLLE_wavespeeds(const Euler3D_Polytropic_pState &Wl,
				    const Euler3D_Polytropic_pState &Wr,
				    const Vector3D &norm_dir);

    //! Returns rotated primitive state aligned from norm_dir to x-direction
    Euler3D_Polytropic_pState Rotate(const Vector3D &norm_dir) const;


    //! Roe Average given 2 primitive states
    static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr);
    //! Roe Average given 2 conservative states
    static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur);

    //! Returns negative waves speeds (eigenvalues) using Harten entropy fix
    static Euler3D_Polytropic_pState lambda_minus(const Euler3D_Polytropic_pState  &lambda_a,
                                                  const Euler3D_Polytropic_pState  &lambda_l,
                                                  const Euler3D_Polytropic_pState  &lambda_r);

    //! Returns positive waves speeds (eigenvalues) using Harten entropy fix
    static Euler3D_Polytropic_pState lambda_plus(const Euler3D_Polytropic_pState  &lambda_a,
                                                 const Euler3D_Polytropic_pState  &lambda_l,
                                                 const Euler3D_Polytropic_pState  &lambda_r);
        
/** @name Boundary Conditions */
/*        ------------------- */
//@{
    //! Return reflected solution state after application of reflection BC
    static Euler3D_Polytropic_pState Reflect(const Euler3D_Polytropic_pState &W, 
                                             const Vector3D &norm_dir);

    //! Return wall solution state after application of moving wall BC
    static Euler3D_Polytropic_pState MovingWall(const Euler3D_Polytropic_pState &Win,
                                                const Euler3D_Polytropic_pState &Wout,
                                                const Vector3D &norm_dir, 
                                                const Vector3D &wall_velocity,
                                                const Vector3D &pressure_gradient,
                                                const int &TEMPERATURE_BC_FLAG);
 
    //! Return wall solution state after application of no-slip BC
    static Euler3D_Polytropic_pState NoSlip(const Euler3D_Polytropic_pState &Win,
                                            const Euler3D_Polytropic_pState &Wout,
                                            const Vector3D &norm_dir,
                                            const Vector3D &pressure_gradient,
                                            const int &TEMPERATURE_BC_FLAG);
//@}
};

/*! 
 * Class: Euler3D_Polytropic_cState
 *
 * \brief Conserved state solution class for 3D Euler equations
 *        governing polytropic gaseous flows.
 *
 *	This class contains the primitive variables and member functions.
 *	It also contains the approximate riemann solvers.
 *
 * Member functions
 *  - rho           -- Return density
 *  - rhov          -- Return mass flux   
 *  - E             -- Return total energy
 *  - g             -- Return specific heat ratio.
 *  - gm1           -- Return g-1
 *  - gm1i          -- Return 1/(g-1).
 *  - R             -- Return gas constant.
 *  - setgas        -- Set gas constants.
 *  - T             -- Return temperature.
 *  - e             -- Return specific internal energy.
 *  - E             -- Return total energy.
 *  - h             -- Return specific enthalpy.
 *  - H             -- Return total enthalpy.
 *  - a             -- Return sound speed.
 *  - a2            -- Return sound speed square.
 *  - M             -- Return Mach number.
 *  - s             -- Return specific entropy.
 *
 * Member operators \n
 *      U -- a conserved solution state \n
 *      c -- a scalar (double)
 *
 *  - U = U;
 *  - c = U[i];
 *  - U = U + U;
 *  - U = U - U;
 *  - c = U * U; (inner product)
 *  - U = c * U;
 *  - U = U * c;
 *  - U = U / c;
 *  - U = U ^ U;
 *  - U = -U;
 *  - U += U;
 *  - U -= U;
 *  - U *= U;
 *  - U /= U;
 *  - U == U;
 *  - U != U;
 *  - cout << U; (output function)
 *  - cin  >> U; (input function)
 */
class Euler3D_Polytropic_cState{
public:
    static int     num_vars;      //!< Number of variables.
    double              rho;      //!< Density.
    Vector3D           rhov;      //!< Momentum.
    double                E;      //!< Total Energy.
    static double         g;      //!< Specific heat ratio.
    static double       gm1;      //!< g-1
    static double      gm1i;      //!< 1/(g-1)
    static double         R;      //!< Gas constant.
    static char*   gas_type;      //!< Gas type

/** @name Constructors and desctructors */
/*        ----------------------------- */
//@{
    //! Creation constructor
    Euler3D_Polytropic_cState(void):
      rho(DENSITY_STDATM) {
      rhov.zero(); E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
    }

    //! Copy constructor
    Euler3D_Polytropic_cState(const Euler3D_Polytropic_cState &U) {
      rho = U.rho; rhov = U.rhov; E = U.E;
    }
    //! Assignment constructor
    Euler3D_Polytropic_cState(const double &d, 
                              const Vector3D &dv, 
                              const double &Etotal):
      rho(d), rhov(dv), E(Etotal) {
    }

    //! Assignment constructor
    Euler3D_Polytropic_cState(const double &d, 
                              const double &dvx, 
                              const double &dvy, 
                              const double &dvz, 
                              const double &Etotal):
      rho(d), E(Etotal) {
      rhov.x = dvx; rhov.y = dvy; rhov.z = dvz;
    }

    //! Assignment constructor
    Euler3D_Polytropic_cState(const Euler3D_Polytropic_pState &W1) {
      rho = W1.rho; rhov = W1.rhov(); E = W1.E();
}
//@}

/** @name Useful operators */
/*        ---------------- */
//@{
    //! Return the number of variables.
    int NumVar(void) {
        return NUM_VAR_EULER3D;
    }

    //! Copy operator.
    void Copy(const Euler3D_Polytropic_cState &U) {
        rho = U.rho;    rhov = U.rhov;      E = U.E;
    }

    //! Vacuum operator.
    void Vacuum(void) {
        rho = ZERO;     rhov = Vector3D_ZERO;   E = ZERO;
    }

    //! Standard atmosphere operator.
    void Standard_Atmosphere(void) {
        rho = DENSITY_STDATM;   rhov.zero();    E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
    }

    //! Check for unphysical state properties.
    int Unphysical_Properties(void) const {
        if (rho <= ZERO || E <= ZERO || e() <= ZERO)
            return 1;
        return 0;
    }
//@}

/** @name Set static variables */ 
/*        -------------------- */
//@{
    void setgas(void);                  //!< Set gas to air 
    void setgas(char *string_ptr);      //!< Set gas constants
//@}
    
/** @name State functions */
/*        --------------- */
//@{
    Vector3D v(void);                   //!< Flow velocity
    double v(const Vector3D &n);        //!< Flow velocity in direction n
    double p(void);                     //!< Pressure
    double T(void);                     //!< Temperature
    double e(void);                     //!< Specific internal energy
    double h(void);                     //!< Specific enthalpy
    double H(void);                     //!< Total enthalpy
    double a(void);                     //!< Sound speed
    double a2(void);                    //!< Sound speed squared
    double M(void);                     //!< Mach number
    double s(void);                     //!< Specific entropy
    double To(void);                    //!< Stagnation temperature
    double po(void);                    //!< Stagnation pressure
    double ao(void);                    //!< Stagnation sound speed
    double ho(void);                    //!< Stagnation enthalpy
//@}

    //same functions with const
    Vector3D v(void) const;
    double v(const Vector3D &n) const;
    double p(void) const;
    double T(void) const;
    double e(void) const;
    double h(void) const;
    double H(void) const;
    double a(void) const;
    double a2(void) const;
    double M(void) const;
    double s(void) const;
    double To(void) const;
    double po(void) const;
    double ao(void) const;
    double ho(void) const;

/** @name Primitive solution state */ 
/*        ------------------------ */
//@{
    //! Convert primitive solution state to conservative solution state
    Euler3D_Polytropic_pState W(void);
    static Euler3D_Polytropic_pState W(const Euler3D_Polytropic_cState &U);
//@}

    //same functions with const
    Euler3D_Polytropic_pState W(void) const;

/** @name Fluxes */
/*        ------ */
//@{
    Euler3D_Polytropic_cState F(void);				//!< x-direction flux
    Euler3D_Polytropic_cState Fx(void);				//!< x-direction flux
    Euler3D_Polytropic_cState Fy(void);				//!< y-direction flux
    Euler3D_Polytropic_cState Fz(void);				//!< z-direction flux
//@}

    //same functions with const
    Euler3D_Polytropic_cState F(void) const;
    static Euler3D_Polytropic_cState F(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState Fx(void) const;
    static Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState Fy(void) const;
    static Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState Fz(void) const;
    static Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_cState &U);

/** @name Flux Jacobians */
/*        -------------- */
//@{
    void dFxdU(DenseMatrix &dFxdU);			//!< x-direction flux Jacobian
    void dFydU(DenseMatrix &dFydU);			//!< y-direction flux Jacobian
    void dFzdU(DenseMatrix &dFzdU);			//!< z-direction flux Jacobian
//@}

    //same functions with const
    void dFxdU(DenseMatrix &dFxdU) const;
    static void dFxdU(DenseMatrix &dFxdU, const Euler3D_Polytropic_cState &U);
    void dFydU(DenseMatrix &dFydU) const;
    static void dFydU(DenseMatrix &dFydU, const Euler3D_Polytropic_cState &U);
    void dFzdU(DenseMatrix &dFzdU) const;
    static void dFzdU(DenseMatrix &dFzdU, const Euler3D_Polytropic_cState &U);

/** @name Solution variable Jacobians. */
/*        ---------------------------- */
//@{
    void dUdW(DenseMatrix &dUdW);		//!< dU/dW
    void dWdU(DenseMatrix &dWdU);		//!< dW/dU
//@}

    // same functions with const
    void dUdW(DenseMatrix &dUdW) const;
    static void dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_cState &U);
    void dWdU(DenseMatrix &dWdU) const;
    static void dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_cState &U);

/** @name Operators. */
/*        ---------- */
//@{
// Index operator. 
    //!  c = U[i]
    double &operator[](int index);		
    const double &operator[](int index) const;

// Binary arithmetic operators. 
    //! U = U + U
    friend Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    //! U = U - U
    friend Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    //! c = U * U (inner product)
    friend double operator					   *(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    //! U = U * c
    friend Euler3D_Polytropic_cState operator *(const Euler3D_Polytropic_cState &U, const double &a);
    //! U = c * U
    friend Euler3D_Polytropic_cState operator *(const double &a, const Euler3D_Polytropic_cState &U);
    //! U = U / c
    friend Euler3D_Polytropic_cState operator /(const Euler3D_Polytropic_cState &U, const double &a);
    //! U = U ^ U
    friend Euler3D_Polytropic_cState operator ^(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);

// Unary arithmetic operators.
    //! U = -U
    friend Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U);

// Shortcut arithmetic operators.
    //! U += U
    Euler3D_Polytropic_cState& operator +=(const Euler3D_Polytropic_cState &U);
    //! U -= U
    Euler3D_Polytropic_cState& operator -=(const Euler3D_Polytropic_cState &U);
    //! U *= U
    Euler3D_Polytropic_cState& operator *=(const double &a);
    //! U /= U
    Euler3D_Polytropic_cState& operator /=(const double &a);

// Relational operators.
    //! U == U
    friend int operator ==(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);
    //! U != U
    friend int operator !=(const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2);

// Input-output operators.
    //! ostream << U
    friend ostream& operator << (ostream &out_file, const Euler3D_Polytropic_cState &U);
    //! istream >> U
    friend istream& operator >> (istream &in_file,  Euler3D_Polytropic_cState &U);
//@}

   //! Returns rotated conserved state aligned with norm_dir
    Euler3D_Polytropic_cState RotateBack(const Vector3D &norm_dir) const;

};

/***************************************************************************************
 * Euler3D_Polytropic_pState member functions                                          *
 ***************************************************************************************/

/***************************************************************************************
 * Euler3D_Polytropic_pState::Euler3D_Polytropic_pState -- Constructor.                *
 ***************************************************************************************/
/*! Creates primitive state from a given conservative state */ 
inline Euler3D_Polytropic_pState::Euler3D_Polytropic_pState(const Euler3D_Polytropic_cState &U1) {
    rho = U1.rho; v = U1.v(); p = U1.p();
}

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Index operators.                                       *
 ***************************************************************************************/
inline double& Euler3D_Polytropic_pState::operator[](int index) {
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

inline const double& Euler3D_Polytropic_pState::operator[](int index) const {
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

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Binary arithmetic operators.                           *
 ***************************************************************************************/
inline Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W1, 
                                            const Euler3D_Polytropic_pState &W2) {
    return (Euler3D_Polytropic_pState(W1.rho+W2.rho,W1.v+W2.v,W1.p+W2.p));
}

inline Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W1, 
                                            const Euler3D_Polytropic_pState &W2) {
    return (Euler3D_Polytropic_pState(W1.rho-W2.rho,W1.v-W2.v, W1.p-W2.p));
}

inline double operator *(const Euler3D_Polytropic_pState &W1, 
                       const Euler3D_Polytropic_pState &W2) {
    return (W1.rho*W2.rho + W1.v.x*W2.v.x + W1.v.y*W2.v.y + W1.v.z*W2.v.z + W1.p*W2.p);
}	

inline Euler3D_Polytropic_pState operator *(const Euler3D_Polytropic_pState &W, 
                                            const double &a) {
    return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z, a*W.p));
}

inline Euler3D_Polytropic_pState operator *(const double &a, 
                                            const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_pState(a*W.rho,a*W.v.x,a*W.v.y, a*W.v.z,a*W.p));
}

inline Euler3D_Polytropic_pState operator /(const Euler3D_Polytropic_pState &W, 
                                            const double &a) {
    return (Euler3D_Polytropic_pState(W.rho/a, W.v.x/a, W.v.y/a, W.v.z/a,W.p/a));
}

inline Euler3D_Polytropic_pState operator ^(const Euler3D_Polytropic_pState &W1, 
                                            const Euler3D_Polytropic_pState &W2) {
    return (Euler3D_Polytropic_pState(W1.rho*W2.rho,
                                      W1.v.x*W2.v.x, 
                                      W1.v.y*W2.v.y, 
                                      W1.v.z*W2.v.z, 
                                      W1.p*W2.p));
}

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Shortcut arithmetic operators.                         *
 ***************************************************************************************/
inline Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator +=(const Euler3D_Polytropic_pState &W) {
    rho += W.rho; 
    v.x += W.v.x; v.y += W.v.y;	v.z += W.v.z;
    p += W.p;
    return *this;
}

inline Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator -=(const Euler3D_Polytropic_pState &W) {
    rho -= W.rho; 
    v.x -= W.v.x; v.y -= W.v.y;	v.z -= W.v.z;
    p -= W.p;
    return *this;
}

inline Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator *=(const double &a) {
    rho *= a; 
    v.x *= a; v.y *= a;	v.z *= a; 
    p *= a;
    return *this;
}

inline Euler3D_Polytropic_pState& Euler3D_Polytropic_pState::operator /=(const double &a) {
    rho /= a; 
    v.x /= a; v.y /= a; v.z /= a;
    p /= a;
    return *this;
}

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Unary arithmetic operators.                            *
 ***************************************************************************************/
inline Euler3D_Polytropic_pState operator -(const Euler3D_Polytropic_pState &W) {
    return (Euler3D_Polytropic_pState(-W.rho,-W.v.x, -W.v.y, -W.v.z, -W.p));
}

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Relational operators.                                  *
 ***************************************************************************************/
inline int operator ==(const Euler3D_Polytropic_pState &W1, 
                       const Euler3D_Polytropic_pState &W2) {
    return (W1.rho == W2.rho && W1.v.x == W2.v.x && W1.v.y == W2.v.y && W1.v.z == W2.v.z && W1.p == W2.p);
}

inline int operator !=(const Euler3D_Polytropic_pState &W1, 
                       const Euler3D_Polytropic_pState &W2) {
    return (W1.rho != W2.rho || W1.v.x != W2.v.x || W1.v.y != W2.v.y || W1.v.z != W2.v.z || W1.p != W2.p);
}

/***************************************************************************************
 * Euler3D_Polytropic_pState -- Input-output operators.                                *
 ***************************************************************************************/
inline ostream& operator << (ostream &out_file, 
                             const Euler3D_Polytropic_pState &W) {
    out_file.setf(ios::scientific);
    out_file << " " << W.rho<< " " << W.v.x << " " << W.v.y 
        << " " << W.v.z << " " << W.p;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream& operator >> (istream &in_file,  
                             Euler3D_Polytropic_pState &W) {
    in_file.setf(ios::skipws);
    in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z 
        >> W.p;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/***************************************************************************************
 * Euler3D_Polytropic_cState member functions                                          *
 ***************************************************************************************/

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Index operators.                                       *
 ***************************************************************************************/
inline double& Euler3D_Polytropic_cState::operator[](int index) {
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

inline const double& Euler3D_Polytropic_cState::operator[](int index) const {
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

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Binary arithmetic operators.                           *
 ***************************************************************************************/
inline Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U1, 
                                            const Euler3D_Polytropic_cState &U2) {
    return (Euler3D_Polytropic_cState(U1.rho+U2.rho,U1.rhov+U2.rhov,U1.E+U2.E));
}

inline Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U1, 
                                            const Euler3D_Polytropic_cState &U2) {
    return (Euler3D_Polytropic_cState(U1.rho-U2.rho,U1.rhov-U2.rhov,U1.E-U2.E));
}

inline double operator *(const Euler3D_Polytropic_cState &U1, 
                         const Euler3D_Polytropic_cState &U2) {
    return (U1.rho*U2.rho+U1.rhov*U2.rhov+U1.E*U2.E);
}

inline Euler3D_Polytropic_cState operator *(const Euler3D_Polytropic_cState &U, 
                                            const double &a) {
    return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
}

inline Euler3D_Polytropic_cState operator *(const double &a, 
                                            const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(a*U.rho,a*U.rhov,a*U.E));
}

inline Euler3D_Polytropic_cState operator /(const Euler3D_Polytropic_cState &U, 
                                            const double &a) {
    return (Euler3D_Polytropic_cState(U.rho/a,U.rhov/a,U.E/a));
}

inline Euler3D_Polytropic_cState operator ^(const Euler3D_Polytropic_cState &U1, 
                                            const Euler3D_Polytropic_cState &U2) {
    return (Euler3D_Polytropic_cState(U1.rho*U2.rho,
                                      U1.rhov.x*U2.rhov.x,
                                      U1.rhov.y*U2.rhov.y,
                                      U1.rhov.z*U2.rhov.z,
                                      U1.E*U2.E));
}

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Shortcut arithmetic operators.                         *
 ***************************************************************************************/
inline Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator +=(const Euler3D_Polytropic_cState &U) {
    rho += U.rho;
    rhov += U.rhov;
    E += U.E;
    return *this;
}

inline Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator -=(const Euler3D_Polytropic_cState &U) {
    rho -= U.rho;
    rhov -= U.rhov;
    E -= U.E;
    return *this;
}

inline Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator *=(const double &a) {
    rho *= a;
    rhov.x *= a; rhov.y *= a; rhov.z *= a;
    E *= a;
    return *this;
}

inline Euler3D_Polytropic_cState& Euler3D_Polytropic_cState::operator /=(const double &a) {
    rho /= a;
    rhov.x /= a; rhov.y /= a; rhov.z /= a;
    E /= a;
    return *this;
}

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Unary arithmetic operators.                            *
 ***************************************************************************************/
inline Euler3D_Polytropic_cState operator -(const Euler3D_Polytropic_cState &U) {
    return (Euler3D_Polytropic_cState(-U.rho,-U.rhov,-U.E));
}

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Relational operators.                                  *
 ***************************************************************************************/
inline int operator ==(const Euler3D_Polytropic_cState &U1, 
                       const Euler3D_Polytropic_cState &U2) {
    return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E);
}

inline int operator !=(const Euler3D_Polytropic_cState &U1, 
                       const Euler3D_Polytropic_cState &U2) {
    return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E);
}

/***************************************************************************************
 * Euler3D_Polytropic_cState -- Input-output operators.                                *
 ***************************************************************************************/
inline ostream& operator << (ostream &out_file, 
                             const Euler3D_Polytropic_cState &U) {
    out_file.setf(ios::scientific);
    out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " "
        << U.rhov.z << " " << U.E;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream& operator >> (istream &in_file, 
                             Euler3D_Polytropic_cState &U) {
    in_file.setf(ios::skipws);
    in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/* -------------------------------------------------------------------------- *
 *                     Useful 3D Euler state constants.                       *
 * -------------------------------------------------------------------------- */
/** @name Usefule 3D Euler state constants. */
//@{
//! Standard atmosphere primitive state
const Euler3D_Polytropic_pState Euler3D_W_STDATM(DENSITY_STDATM, Vector3D_ZERO, PRESSURE_STDATM);

//! Vacuum primitive state
const Euler3D_Polytropic_pState Euler3D_W_VACUUM(ZERO, Vector3D_ZERO, ZERO);

//! Standard atmosphere conservative state
const Euler3D_Polytropic_cState Euler3D_U_STDATM(Euler3D_W_STDATM);

//! Vacuum conservative state
const Euler3D_Polytropic_cState Euler3D_U_VACUUM(Euler3D_W_VACUUM);
//@}

#endif /* _EULER3D_POLYTROPIC_STATE_INCLUDED  */

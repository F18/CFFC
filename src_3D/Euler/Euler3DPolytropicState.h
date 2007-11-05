/*! \file Euler3DPolytropicState.h
 * 	\brief	This file defines the state classes for 
 *			Euler3D Polytropic flows.
 */

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

/*! Number of variables in Euler 3D*/
#define	NUM_VAR_EULER3D    5


/*! \brief This is the primitive state class for Euler 3D Polytropic flows.
 *
 *	This class contains the primitive variables and memberfunctions.
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
 *  - W = +W;
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

/*
 * Constructors
 * ------------
 */
        
    //! Creation constructor 
    Euler3D_Polytropic_pState();

    //! Copy constructor 
    Euler3D_Polytropic_pState(const Euler3D_Polytropic_pState &W);

    //! Assignment constructor 
    Euler3D_Polytropic_pState(const double &d, 
                              const Vector3D &V, 
                              const double &pre);
    //! Assignment constructor 
    Euler3D_Polytropic_pState(const double &d, 
                              const double &vx, const double &vy, const double &vz, 
                              const double &pre);
    //! Assignment constructor 
    Euler3D_Polytropic_pState(const Euler3D_Polytropic_cState &U);
        
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
        
/** @name State functions */
/*        --------------- */
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
//! Convert primitive solution state to conservative solution state
//@{
    Euler3D_Polytropic_cState U(void);
    static Euler3D_Polytropic_cState U(const Euler3D_Polytropic_pState &W);
//@}
    // same functions with const
    Euler3D_Polytropic_cState U(void) const;

        
/** @name Fluxes */
/*        ------ */
//@{
    Euler3D_Polytropic_cState Fx(void);         //!< x-direction Flux
    Euler3D_Polytropic_cState Fy(void);         //!< y-direction Flux
    Euler3D_Polytropic_cState Fz(void);         //!< z-direction Flux
//@}
    //same functions with const
    Euler3D_Polytropic_cState Fx(void) const;
    static Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_cState Fy(void) const;
    static Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_cState Fz(void) const;
    static Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_pState &W);
        
        
/** @name Flux Jacobians */
/*        -------------- */
    //@{
    void dFxdU(DenseMatrix &dFxdU);         //!< x-direction Flux Jacobian
    void dFydU(DenseMatrix &dFydU);         //!< y-direction Flux Jacobian
    void dFzdU(DenseMatrix &dFzdU);         //!< z-direction Flux Jacobian
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
    void dUdW(DenseMatrix &dUdW);         //!< dUdW
    void dWdU(DenseMatrix &dWdU);         //!< dWdU
    //@}
    // same functions with const
    void dUdW(DenseMatrix &dUdW) const;
    static void dUdW(DenseMatrix &dUdW, const Euler3D_Polytropic_pState &W);
    void dWdU(DenseMatrix &dWdU) const;
    static void dWdU(DenseMatrix &dWdU, const Euler3D_Polytropic_pState &W);
        
/** @name Eigenvalues */
/*        ----------- */
    //@{	
    Euler3D_Polytropic_pState lambda_x(void);   //!< x-direction eigenvalues
    double lambda_x(int index);                 //!< x-direction eigenvalues

    Euler3D_Polytropic_pState lambda_y(void);   //!< y-direction eigenvalues
    double lambda_y(int index);                 //!< y-direction eigenvalues

    Euler3D_Polytropic_pState lambda_z(void);   //!< z-direction eigenvalues
    double lambda_z(int index);                 //!< z-direction eigenvalues
    //@}
    //same functions with const
    Euler3D_Polytropic_pState lambda_x(void) const;
    static Euler3D_Polytropic_pState lambda_x(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState lambda_y(void) const;
    static Euler3D_Polytropic_pState lambda_y(const Euler3D_Polytropic_pState &W);
    Euler3D_Polytropic_pState lambda_z(void) const;
    static Euler3D_Polytropic_pState lambda_z(const Euler3D_Polytropic_pState &W);
    double lambda_x(int index) const;
    double lambda_y(int index) const;
    double lambda_z(int index) const;




/** @name Conserved right eigenvector */
/*        --------------------------- */
    //@{
    Euler3D_Polytropic_cState rc_x(const int &index);   //!< x-direction conserved right eigenvector
    Euler3D_Polytropic_cState rc_y(const int &index);   //!< y-direction conserved right eigenvector
    Euler3D_Polytropic_cState rc_z(const int &index);   //!< z-direction conserved right eigenvector
    //@}
    //same functions with const
    Euler3D_Polytropic_cState rc_x(const int &index) const;
    Euler3D_Polytropic_cState rc_y(const int &index) const;
    Euler3D_Polytropic_cState rc_z(const int &index) const;


/** @name  Primitive left eigenvector */
/*         -------------------------- */
    //@{
    Euler3D_Polytropic_pState lp_x(const int &index);   //!< x-direction primitive left eigenvector
    Euler3D_Polytropic_pState lp_y(const int &index);   //!< y-direction primitive left eigenvector
    Euler3D_Polytropic_pState lp_z(const int &index);   //!< z-direction primitive left eigenvector
    //@}
    //same functions with const
    Euler3D_Polytropic_pState lp_x(const int &index) const;
    Euler3D_Polytropic_pState lp_y(const int &index) const;
    Euler3D_Polytropic_pState lp_z(const int &index) const;

/** @name Primitive right eigenvector */
/*        --------------------------- */
    //@{
    Euler3D_Polytropic_pState rp_x(const int &index);   //!< x-direction primitive right eigenvector
    Euler3D_Polytropic_pState rp_y(const int &index);   //!< y-direction primitive right eigenvector
    Euler3D_Polytropic_pState rp_z(const int &index);   //!< z-direction primitive right eigenvector
    //@}
    //same functions with const
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
    //! W = +W
    friend Euler3D_Polytropic_pState operator +(const Euler3D_Polytropic_pState &W);
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
        
/** @name Approximate Riemann solvers */
/*        --------------------------- */
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

    //! Returns rotated primitive state aligned with norm_dir
    static Euler3D_Polytropic_pState Rotate(const Vector3D &norm_dir);

    //! Roe Average given 2 primitive states
    static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_pState &Wl,
                                                const Euler3D_Polytropic_pState &Wr);
    //! Roe Average given 2 conservative states
    static Euler3D_Polytropic_pState RoeAverage(const Euler3D_Polytropic_cState &Ul,
                                                const Euler3D_Polytropic_cState &Ur);

    friend Euler3D_Polytropic_pState HartenFixNeg(const Euler3D_Polytropic_pState  &lambda_a,
                                                  const Euler3D_Polytropic_pState  &lambda_l,
                                                  const Euler3D_Polytropic_pState  &lambda_r);

    friend Euler3D_Polytropic_pState HartenFixPos(const Euler3D_Polytropic_pState  &lambda_a,
                                                  const Euler3D_Polytropic_pState  &lambda_l,
                                                  const Euler3D_Polytropic_pState  &lambda_r);
        
/** @name Boundary Conditions */
/*        ------------------- */ 
    //@{
    static Euler3D_Polytropic_pState Reflect(const Euler3D_Polytropic_pState &W, 
                                             const Vector3D &norm_dir);
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
    //@}
        
};


/*! \brief This is the conserved state class for Euler 3D Polytropic flows.
 *
 * This class contains the conserved variables and memberfunctions
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


/*
 * Constructors
 * ------------
 */

    //! Creation constructor
    Euler3D_Polytropic_cState();

    //! Copy constructor
    Euler3D_Polytropic_cState(const Euler3D_Polytropic_cState &U);

    //! Assignment constructor
    Euler3D_Polytropic_cState(const double &rho, 
                              const Vector3D &rhov, 
                              const double &E);

    //! Assignment constructor
    Euler3D_Polytropic_cState(const double &rho, 
                              const double &rhovx, const double &rhovy, const double &rhovz, 
                              const double &E);

    //! Assignment constructor
    Euler3D_Polytropic_cState(const Euler3D_Polytropic_pState &W);

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
    Euler3D_Polytropic_cState Fx(void);				//!< x-direction Flux
    Euler3D_Polytropic_cState Fy(void);				//!< y-direction Flux
    Euler3D_Polytropic_cState Fz(void);				//!< z-direction Flux
    //@}
    //same functions with const
    Euler3D_Polytropic_cState Fx(void) const;
    static Euler3D_Polytropic_cState Fx(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState Fy(void) const;
    static Euler3D_Polytropic_cState Fy(const Euler3D_Polytropic_cState &U);
    Euler3D_Polytropic_cState Fz(void) const;
    static Euler3D_Polytropic_cState Fz(const Euler3D_Polytropic_cState &U);


/** @name Flux Jacobians */
/*        -------------- */
    //@{
    void dFxdU(DenseMatrix &dFxdU);			//!< x-direction Flux Jacobian
    void dFydU(DenseMatrix &dFydU);			//!< y-direction Flux Jacobian
    void dFzdU(DenseMatrix &dFzdU);			//!< z-direction Flux Jacobian
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
    void dUdW(DenseMatrix &dUdW);		//!< dUdW
    void dWdU(DenseMatrix &dWdU);		//!< dWdU
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
    //! U = +U
    friend Euler3D_Polytropic_cState operator +(const Euler3D_Polytropic_cState &U);
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
    
};




/* -------------------------------------------------------------------------- *
 *                     Useful 3D Euler state constants.                       *
 * -------------------------------------------------------------------------- */
/** @name Usefule 3D Euler state constants. */
//@{
//! Standard atmosphere primitive state
const Euler3D_Polytropic_pState Euler3D_W_STDATM(DENSITY_STDATM,
				      Vector3D_ZERO, PRESSURE_STDATM);
//! Vacuum primitive state
const Euler3D_Polytropic_pState Euler3D_W_VACUUM(ZERO, Vector3D_ZERO, ZERO);
//! Standard atmosphere conservative state
const Euler3D_Polytropic_cState Euler3D_U_STDATM(Euler3D_W_STDATM);
//! Vacuum conservative state
const Euler3D_Polytropic_cState Euler3D_U_VACUUM(Euler3D_W_VACUUM);
//@}

#endif /* _EULER3D_POLYTROPIC_STATE_INCLUDED  */

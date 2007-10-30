/* Gaussian2DState.h:  Header file defining 2D Gaussian Solution State Classes. */

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#define _GAUSSIAN2D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, 2D vector, and gas constant header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

#ifndef _THIRD_ORDER_TENSOR2D_INCLUDED
#include "../Math/Third_order_tensor2D.h"
#endif //_THIRD_ORDER_TENSOR2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _SOLID_CONSTANTS_INCLUDED
#include "../Physics/SolidConstants.h"
#endif // _SOLID_CONSTANTS_INCLUDED //needed for embeddedboundaries2D

/* Define the classes. */

#define	NUM_VAR_GAUSSIAN2D    8
#define GAUSSIAN_MONATOMIC    1
#define GAUSSIAN_DIATOMIC     2

#define _GAUSSIAN_HEAT_TRANSFER_


class Gaussian2D_cState;
class Gaussian2D_Input_Parameters;

/********************************************************
 * Class: Gaussian2D_pState                             *
 *                                                      *
 * Member functions                                     *
 *     d        -- Return density.                      *
 *     v        -- Return flow velocity.                *
 *     p        -- Return pressure.                     *
 *     q        -- third-order heat flux tensor         *
 *                 (for regularized 10-moment)          *
 *     atoms    -- Return number of atoms in a molecule *
 *     M        -- Return Molar mass                    *
 *     setgas   -- Set gas constants.                   *
 *     T        -- Return temperature.                  *
 *     axx      -- Return pxx/rho.                      *
 *     ayy      -- Return pyy/rho.                      *
 *     sound    -- Return sound speed.                  *
 *     R        -- Return specific gas constant.        *
 *     Gamma    -- Return gamma for gas.                *
 *     Cv       -- Return specific heat at constant V.  *
 *     Cp       -- Return specific heat at constant P.  *
 *     K        -- Return thermal conductivity.         *
 *     dv       -- Return momentum.                     *
 *     U        -- Return conserved solution state.     *
 *     F        -- Return x-direction solution flux.    *
 *     Fx       -- Return x-direction solution flux.    *
 *     Fy       -- Return y-direction solution flux.    *
 *     Fn       -- Return n-direction solution flux.    *
 *     Gx       -- Return x-direction elliptic flux.    *
 *     Gy       -- Return y-direction elliptic flux.    *
 *     lambda   -- Return x-direction eigenvalue(s).    *
 *     lambda_x -- Return x-direction eigenvalue(s).    *
 *     lambda_y -- Return y-direction eigenvalue(s).    *
 *     rp       -- Return primitive right eigenvector   *
 *                 (x-direction).                       *
 *     rp_x     -- Return primitive right eigenvector   *
 *                 (x-direction).                       *
 *     rp_y     -- Return primitive right eigenvector   *
 *                 (y-direction).                       *
 *     rc       -- Return conserved right eigenvector   *
 *                 (x-direction).                       *
 *     rc_x     -- Return conserved right eigenvector   *
 *                 (x-direction).                       *
 *     rc_y     -- Return conserved right eigenvector   *
 *                 (y-direction).                       *
 *     lp       -- Return primitive left eigenvector    *
 *                 (x-direction).                       *
 *     lp_x     -- Return primitive left eigenvector    *
 *                 (x-direction).                       *
 *     lp_y     -- Return primitive left eigenvector    *
 *                 (y-direction).                       *
 *     S        -- Return axisymmetric source term      *
 *                 vector.                              *
 *     tt       -- translational relaxation time        *
 *     tr       -- rotational relaxation time           *
 *     gt       -- temperature jump distance            *
 *                 for temperature slip BC              *
 *     relax    -- relaxes state towards thermodynamic  *
 *                 equilibrium (BGK source terms)       *
 *     set_temperature_d                                *
 *              -- change density so that the state has *
 *                 the desired temperature              *
 *     set_state_from_ips                               *
 *              -- set ics relating to temperature and  *
 *                 mach number etc.                     *
 *                                                      *
 * Member operators                                     *
 *      W -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * W = W;                                               *
 * c = W[i];                                            *
 * W = W + W;                                           *
 * W = W - W;                                           *
 * c = W * W; (inner product)                           *
 * W = c * W;                                           *
 * W = W * c;                                           *
 * W = W / c;                                           *
 * W = W ^ W; (my useful product)                       *
 * W = +W;                                              *
 * W = -W;                                              *
 * W += W;                                              *
 * W -= W;                                              *
 * W == W;                                              *
 * W != W;                                              *
 * cout << W; (output function)                         *
 * cin  >> W; (input function)                          *
 *                                                      *
 ********************************************************/
class Gaussian2D_pState{
  private:
  public:
    double                    d;   // Density.
    Vector2D                  v;   // Flow velocity (2D vector).
    Tensor2D                  p;   // Pressure tensor.
    double                 erot;   // rotational energy
#ifdef _GAUSSIAN_HEAT_TRANSFER_
    Third_order_tensor2D      q;   // third-order heat tensor
#endif
    static double             M;   // Molar Weight
    static int            atoms;   // Monatomic or Diatomic
    static int              gas;   // Gas type
    static double       alpha_m;   // Momentum accommodation coefficient for gas/solid boundary
    static double       alpha_t;   // Thermal accommodation coefficient for gas/solid boundary
    static double         omega;   // Viscosity exponent
    static double        mu_not;   // Reference viscosity
    static double            pr;   // Prandtl number
                                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_pState(void) {
       d = DENSITY_STDATM; v.zero(); p.xx = PRESSURE_STDATM;
       p.xy = ZERO; p.yy = PRESSURE_STDATM; p.zz = PRESSURE_STDATM;
       if(atoms==2) {erot = PRESSURE_STDATM;}
       else{erot = 0.0;}
    }

    Gaussian2D_pState(const Gaussian2D_pState &W) {
       d = W.d; v = W.v; p = W.p; erot = W.erot;
    }

    Gaussian2D_pState(const Gaussian2D_cState &U);

    Gaussian2D_pState(const double &rho,
	           const Vector2D &V,
	           const double &pre) {
       d = rho; v = V; p.xx = pre;
       p.xy = ZERO; p.yy = pre; p.zz = pre;
       if(atoms==2) {erot = pre;}
       else{erot = 0.0;}
    }

    Gaussian2D_pState(const double &rho,
	           const double &vx,
	           const double &vy,
	           const double &pre) {
       d = rho; v.x = vx; v.y = vy; p.xx = pre;
       p.xy = ZERO; p.yy = pre; p.zz = pre;
       if(atoms==2) {erot = pre;}
       else{erot = 0.0;}
    }

    Gaussian2D_pState(const double &rho,
		      const Vector2D &V,
		      const Tensor2D &pre,
		      const double &energyrot) {
      d = rho; v = V; p = pre; erot = energyrot;
    }

    Gaussian2D_pState(const double &rho,
		      const double &vx,
		      const double &vy,
		      const Tensor2D &pre,
		      const double &energyrot) {
      d = rho; v.x = vx; v.y = vy; p = pre; 
      erot = energyrot;
    }


    Gaussian2D_pState(const double &rho,
	           const double &vx,
	           const double &vy,
	           const double &pxx,
                   const double &pxy,
                   const double &pyy,
                   const double &pzz,
                   const double &energyrot) {
       d = rho; v.x = vx; v.y = vy; p.xx = pxx;
       p.xy = pxy; p.yy = pyy; p.zz = pzz;
       erot = energyrot;
    }

    
    /* Destructor. */
    // ~Gaussian2D_pState(void);
    // Use automatically generated destructor.

    /* Return the number of variables. */
    int NumVar(void) { return NUM_VAR_GAUSSIAN2D; }

    /* Vacuum operator. */
    void Vacuum(void) {
       d = ZERO; v = Vector2D_ZERO; p.xx = ZERO;
       p.xy = ZERO; p.yy = ZERO; p.zz = ZERO; erot = ZERO;
#ifdef _GAUSSIAN_HEAT_TRANSFER_
       q.zero();
#endif
    }

    /* Standard atmosphere operator. */
    void Standard_Atmosphere(void) {
       d = DENSITY_STDATM; v.zero(); p.xx = PRESSURE_STDATM;
       p.xy = ZERO; p.yy = PRESSURE_STDATM; p.zz = PRESSURE_STDATM; erot = PRESSURE_STDATM;
#ifdef _GAUSSIAN_HEAT_TRANSFER_
       q.zero();
#endif
    }


    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Thermodynamic Pressure */

    double pressure(void);
    double pressure(void) const;

    /* Determinant of Pressure tensor */

    double DetP(void);
    double DetP(void) const;

    /* Validity check and adjustment */
    int invalid() const;
    int Unphysical_Properties(void) const {return invalid();}
    void make_valid();

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Total energy. */
    Tensor2D E(void);
    Tensor2D E(void) const;

    /* Specific enthalpy. */
    Tensor2D h(void);
    Tensor2D h(void) const;

    /* Sound speeds. */
    double axx(void);
    double axx(void) const;

    double ayy(void);
    double ayy(void) const;

    double sound(void);
    double sound(void) const;

    /* Other gas properties */
    double R(void) const;
    double Gamma(void) const;
    double Cv(void) const;
    double Cp(void) const;
    double K(void) const;

    /* Momentum. */
    Vector2D dv(void);
    Vector2D dv(void) const;
    double dv(const Vector2D &n);
    double dv(const Vector2D &n) const;

    /* Visocsities. */
    double viscosity(void);
    double viscosity(void) const;
    double mu(void) {return viscosity();}
    double mu(void) const {return viscosity();}
    double nu(void);
    double nu(void) const;
    double bulk_viscosity(void);
    double bulk_viscosity(void) const;

    /* Burning rate. */
    double burningrate(void) const; //needed for compatibility with embeddedboundaries2D

    /* Mean Free Path */
    double mfp(void);

    /* Conserved solution state. */
    Gaussian2D_cState U(void);
    Gaussian2D_cState U(void) const;
    Gaussian2D_cState U(const Gaussian2D_pState &W);
    friend Gaussian2D_cState U(const Gaussian2D_pState &W);
    
    /* Solution flux (x-direction). */
    Gaussian2D_cState F(void);
    Gaussian2D_cState F(void) const;
    Gaussian2D_cState F(const Gaussian2D_pState &W);
    friend Gaussian2D_cState F(const Gaussian2D_pState &W);
    Gaussian2D_cState F(const Vector2D &V) const;

    Gaussian2D_cState Fx(void);
    Gaussian2D_cState Fx(void) const;
    Gaussian2D_cState Fx(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fx(const Gaussian2D_pState &W);

    /* Solution flux (y-direction). */
    Gaussian2D_cState Fy(void);
    Gaussian2D_cState Fy(void) const;
    Gaussian2D_cState Fy(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fy(const Gaussian2D_pState &W);

    /* Solution flux (n-direction). */
    Gaussian2D_cState Fn(void);
    Gaussian2D_cState Fn(void) const;
    Gaussian2D_cState Fn(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fn(const Gaussian2D_pState &W);

    // Elliptic Heat terms
    Gaussian2D_cState Gx(const Gaussian2D_pState &dWdx) const;
    Gaussian2D_cState Gy(const Gaussian2D_pState &dWdy) const;

    // Compute Heat terms
    void ComputeHeatTerms(const Gaussian2D_pState &dWdx,
			  const Gaussian2D_pState &dWdy,
			  const Vector2D &X,
			  const int &Axisymmetric);

    /* Eigenvalue(s) (x-direction). */
    Gaussian2D_pState lambda(void);
    Gaussian2D_pState lambda(void) const;
    Gaussian2D_pState lambda(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda(const Gaussian2D_pState &W);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const Gaussian2D_pState &W, int index);

    Gaussian2D_pState lambda_x(void);
    Gaussian2D_pState lambda_x(void) const;
    Gaussian2D_pState lambda_x(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda_x(const Gaussian2D_pState &W);
    double lambda_x(int index);
    double lambda_x(int index) const;
    friend double lambda_x(const Gaussian2D_pState &W, int index);
    Gaussian2D_pState lambda_x(const Vector2D &V) const;

    /* Eigenvalue(s) (y-direction). */
    Gaussian2D_pState lambda_y(void);
    Gaussian2D_pState lambda_y(void) const;
    Gaussian2D_pState lambda_y(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda_y(const Gaussian2D_pState &W);
    double lambda_y(int index);
    double lambda_y(int index) const;
    friend double lambda_y(const Gaussian2D_pState &W, int index);

    /* Conserved right eigenvector (x-direction). */
    Gaussian2D_cState rc(int index);
    Gaussian2D_cState rc(int index) const;
    friend Gaussian2D_cState rc(const Gaussian2D_pState &W, int index);

    Gaussian2D_cState rc_x(int index);
    Gaussian2D_cState rc_x(int index) const;
    friend Gaussian2D_cState rc_x(const Gaussian2D_pState &W, int index);

   /* Conserved right eigenvector (y-direction). */
    Gaussian2D_cState rc_y(int index);
    Gaussian2D_cState rc_y(int index) const;
    friend Gaussian2D_cState rc_y(const Gaussian2D_pState &W, int index);

    /* Primitive left eigenvector (x-direction). */
    Gaussian2D_pState lp(int index);
    Gaussian2D_pState lp(int index) const;
    friend Gaussian2D_pState lp(const Gaussian2D_pState &W, int index);

    Gaussian2D_pState lp_x(int index);
    Gaussian2D_pState lp_x(int index) const;
    friend Gaussian2D_pState lp_x(const Gaussian2D_pState &W, int index);

    /* Primitive left eigenvector (y-direction). */
    Gaussian2D_pState lp_y(int index);
    Gaussian2D_pState lp_y(int index) const;
    friend Gaussian2D_pState lp_y(const Gaussian2D_pState &W, int index);

    /* Source terms */
    void relax(double deltat, int stage, const Gaussian2D_pState &W);
    int analytically_inverted_relaxation() { //this is needed for embeddedboundaries with gaussian2D
      return 1;
    }

    //translational and rotational relaxation times
    double tt() const;
    double tr() const;

    //temperature slip distance
    double gt() const;

    //set temperature and set state from ics
    void set_temperature_d(double temperature);
    void set_state_from_ips(Gaussian2D_Input_Parameters &IP); //in Gaussian2D_State.cc

    /* Source vector (axisymmetric terms). */
    Gaussian2D_cState S(const Vector2D &X);
    Gaussian2D_cState S(const Vector2D &X) const;
    friend Gaussian2D_cState S(const Gaussian2D_pState &W, const Vector2D &X);

    void dSdU(DenseMatrix &dSdU,
	      const Vector2D &X,
	      const Gaussian2D_pState &dWdx,
	      const Gaussian2D_pState &dWdy,
	      const int &Axisymmetric) {}

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (p.xx);
        case 5 :
	  return (p.xy);
        case 6 :
	  return (p.yy);
        case 7 :
	  return (p.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (p.xx);
        case 5 :
	  return (p.xy);
        case 6 :
	  return (p.yy);
        case 7 :
	  return (p.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }

    // Assignment operator. I added it due to a bug in the icc compiler
    Gaussian2D_pState& operator =(const Gaussian2D_pState &W);

    /* Binary arithmetic operators. */
    friend Gaussian2D_pState operator +(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState operator -(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend double operator *(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState operator *(const Gaussian2D_pState &W, const double &a);
    friend Gaussian2D_pState operator *(const double &a, const Gaussian2D_pState &W);
    friend Gaussian2D_pState operator /(const Gaussian2D_pState &W, const double &a);
    friend Gaussian2D_pState operator ^(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    
    /* Unary arithmetic operators. */
    friend Gaussian2D_pState operator +(const Gaussian2D_pState &W);
    friend Gaussian2D_pState operator -(const Gaussian2D_pState &W);

    /* Shortcut arithmetic operators. */
    friend Gaussian2D_pState &operator +=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState &operator -=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    Gaussian2D_pState &operator *=(const double &a);
    Gaussian2D_pState &operator /=(const double &a);

    /* Relational operators. */
    friend int operator ==(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend int operator !=(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Gaussian2D_pState &W);
    friend istream &operator >> (istream &in_file,  Gaussian2D_pState &W);

    void output_labels(ostream &out_file) {
      out_file << "\"rho\" \\ \n"
	       << "\"u\" \\ \n"
	       << "\"v\" \\ \n"
	       << "\"pxx\" \\ \n"
	       << "\"pxy\" \\ \n"
	       << "\"pyy\" \\ \n"
	       << "\"pzz\" \\ \n"
	       << "\"erot\" \\ \n"
	       << "\"Pressure\" \\ \n"
	       << "\"Mach\" \\ \n";
    }

    void output_data(ostream &out_file, double dummy1, double dummy2) { //dummies needed for compatibility with NS turbulent
      out_file << " " << d << " " << v.x << " " << v.y << " " << p.xx
	       << " " << p.xy << " " << p.yy << " " << p.zz << " " << erot
	       << " " << pressure() << " " << (sqrt(sqr(v.x)+sqr(v.y))/sound());
    }


    double a(void) const {return 0.0;} //needed for compatibility with embeddedboundaries2D

};

/********************************************************
 * Class: Gaussian2D_cState                             *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      E       -- Return total energy.                 *
 *      q       -- third-order heat flux tensor         *
 *                 (for regularized 10-moment)          *
 *      atoms   -- Return number of atoms in a molecule *
 *      M       -- Return molar mass                    *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      W       -- Return primitive solution state.     *
 *      F       -- Return x-direction solution flux.    *
 *      Fx      -- Return x-direction solution flux.    *
 *      Fy      -- Return y-direction solution flux.    *
 *      Fn      -- Return n-direction solution flux.    *
 *      S       -- Return axisymmetric source term      *
 *                 vector.                              *
 *      Gx      -- Return x-direction elliptic flux.    *
 *      Gy      -- Return y-direction elliptic flux.    *
 *                                                      *
 *     ComputeHeatTerms                                 *
 *              -- compute Heat Transfer terms          *
 *                                                      *
 *                                                      *
 * Member operators                                     *
 *      U -- a conserved solution state                 *
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
 * U = U ^ U; (my useful product)                       *
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
class Gaussian2D_cState{
  private:
  public:
    double                   d;   // Density.
    Vector2D                dv;   // Momentum.
    Tensor2D                 E;   // Total Energy.
#ifdef _GAUSSIAN_HEAT_TRANSFER_
    Third_order_tensor2D     q;   // third-order heat tensor
#endif
    double                erot;   // Rotational Energy
    static double            M;   // Molar Weight
    static int           atoms;   // Monatomic or Diatomic
    static int             gas;   // gas type
    static double      alpha_m;   // Momentum accommodation coefficient for gas/solid boundary
    static double      alpha_t;   // Thermal accommodation coefficient for gas/solid boundary
    static double           pr;   // Prandtl number
	                          // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_cState(void) {
       d = DENSITY_STDATM; dv.zero(); E.xx = PRESSURE_STDATM;
       E.xy = ZERO; E.yy = PRESSURE_STDATM; E.zz = PRESSURE_STDATM;
       if(atoms==2) {erot = PRESSURE_STDATM;}
       else{erot = 0.0;}

    }

    Gaussian2D_cState(const Gaussian2D_cState &U) {
       d = U.d; dv = U.dv; E = U.E; erot = U.erot;
    }

    Gaussian2D_cState(const Gaussian2D_pState &W);

    Gaussian2D_cState(const double &rho,
	           const Vector2D &rhoV,
	           const Tensor2D &Etotal,
                   const double &energyrot) {
       d = rho; dv = rhoV; E = Etotal;
       erot = energyrot;
    }

    Gaussian2D_cState(const double &rho,
	           const double &rhovx,
	           const double &rhovy,
	           const double &exx,
                   const double &exy,
                   const double &eyy,
                   const double &ezz,
                   const double &energyrot) {
       d = rho; dv.x = rhovx; dv.y = rhovy; E.xx = exx;
       E.xy = exy, E.yy = eyy; E.zz = ezz; erot = energyrot;
    }
    
    /* Destructor. */
    // ~Gaussian2D_cState(void);
    // Use automatically generated destructor.

    /* Return the number of variables. */
    int NumVar(void) { return NUM_VAR_GAUSSIAN2D; }

    /* Vacuum operator. */
    void Vacuum(void) {
       d = ZERO; dv = Vector2D_ZERO; E.xx = ZERO;
       E.xy = ZERO; E.yy = ZERO; E.zz = ZERO; erot = ZERO; 
#ifdef _GAUSSIAN_HEAT_TRANSFER_
       q.zero();
#endif
    }

    /* Standard atmosphere operator. */
    void Standard_Atmosphere(void) {
       d = DENSITY_STDATM; dv.zero(); E.xx = PRESSURE_STDATM;
       E.xy = ZERO; E.yy = PRESSURE_STDATM; E.zz = PRESSURE_STDATM; erot = PRESSURE_STDATM;
#ifdef _GAUSSIAN_HEAT_TRANSFER_
       q.zero();
#endif
    }

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Flow velocity. */
    Vector2D v(void);
    Vector2D v(void) const;
    double v(const Vector2D &n);
    double v(const Vector2D &n) const;

    /* Pressure. */
    Tensor2D p(void);
    Tensor2D p(void) const;

    /* Thermodynamic Pressure */
    double pressure(void) const;

    /* Validity check */
    int invalid()const;
    int Unphysical_Properties(void) const {return invalid();}

    /* Burning rate. */
    double burningrate(void) const; //used for embeddedboundaries2D

    /* Primitive solution state. */
    Gaussian2D_pState W(void);
    Gaussian2D_pState W(void) const;
    Gaussian2D_pState W(const Gaussian2D_cState &U);
    friend Gaussian2D_pState W(const Gaussian2D_cState &U);
    
    /* Solution flux (x-direction). */
    Gaussian2D_cState F(void);
    Gaussian2D_cState F(void) const;
    Gaussian2D_cState F(const Gaussian2D_cState &U);
    friend Gaussian2D_cState F(const Gaussian2D_cState &U);
    Gaussian2D_cState F(const Vector2D &V) const;

    Gaussian2D_cState Fx(void);
    Gaussian2D_cState Fx(void) const;
    Gaussian2D_cState Fx(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fx(const Gaussian2D_cState &U);

    /* Solution flux (y-direction). */
    Gaussian2D_cState Fy(void);
    Gaussian2D_cState Fy(void) const;
    Gaussian2D_cState Fy(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fy(const Gaussian2D_cState &U);

    /* Solution flux (n-direction). */
    Gaussian2D_cState Fn(void);
    Gaussian2D_cState Fn(void) const;
    Gaussian2D_cState Fn(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fn(const Gaussian2D_cState &U);

    // Elliptic Heat terms
    Gaussian2D_cState Gx(const Gaussian2D_pState &dWdx) const;
    Gaussian2D_cState Gy(const Gaussian2D_pState &dWdy) const;

    // Compute Heat terms
    void ComputeHeatTerms(const Gaussian2D_pState &dWdx,
			  const Gaussian2D_pState &dWdy,
			  const Vector2D &X,
			  const int &Axisymmetric);

    /* Source vector (axisymmetric terms). */
    Gaussian2D_cState S(const Vector2D &X);
    Gaussian2D_cState S(const Vector2D &X) const;
    friend Gaussian2D_cState S(const Gaussian2D_cState &U, const Vector2D &X);

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (E.xx);
        case 5 :
	  return (E.xy);
        case 6 :
	  return (E.yy);
        case 7 :
	  return (E.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (E.xx);
        case 5 :
	  return (E.xy);
        case 6 :
	  return (E.yy);
        case 7 :
	  return (E.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }

    // Assignment operator.
    Gaussian2D_cState& operator =(const Gaussian2D_cState &U);

    /* Binary arithmetic operators. */
    friend Gaussian2D_cState operator +(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState operator -(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend double operator *(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState operator *(const Gaussian2D_cState &U, const double &a);
    friend Gaussian2D_cState operator *(const double &a, const Gaussian2D_cState &U);
    friend Gaussian2D_cState operator /(const Gaussian2D_cState &U, const double &a);
    friend Gaussian2D_cState operator ^(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);

    /* Unary arithmetic operators. */
    friend Gaussian2D_cState operator +(const Gaussian2D_cState &U);
    friend Gaussian2D_cState operator -(const Gaussian2D_cState &U);

    /* Shortcut arithmetic operators. */
    friend Gaussian2D_cState &operator +=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState &operator -=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    Gaussian2D_cState &operator *=(const double &a);
    Gaussian2D_cState &operator /=(const double &a);

    /* Relational operators. */
    friend int operator ==(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend int operator !=(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Gaussian2D_cState &U);
    friend istream &operator >> (istream &in_file,  Gaussian2D_cState &U);
    
};

/********************************************************
 * Gaussian2D_pState::setgas -- Assign gas constants.   *
 ********************************************************/
inline void Gaussian2D_pState::setgas(void) {
  M = MOLE_WT_AIR;
  atoms = GAUSSIAN_DIATOMIC;
  gas = GAS_AIR;
  omega = OMEGA_AIR;
  mu_not = MU_NOT_AIR;
}

inline void Gaussian2D_pState::setgas(char *string_ptr) {
  if (strcmp(string_ptr, "AIR") == 0) {
    M = MOLE_WT_AIR;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_AIR;
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
  } else if (strcmp(string_ptr, "A") == 0) {
    M = MOLE_WT_A;
    atoms = GAUSSIAN_MONATOMIC;
    gas = GAS_A;
    omega = OMEGA_A;
    mu_not = MU_NOT_A;
  } else if (strcmp(string_ptr, "CO") == 0) {
    M = MOLE_WT_CO;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_CO;
    omega = OMEGA_CO;
    mu_not = MU_NOT_CO;
  } else if (strcmp(string_ptr, "H2") == 0) {
    M = MOLE_WT_H2;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_H2;
    omega = OMEGA_H2;
    mu_not = MU_NOT_H2;
  } else if (strcmp(string_ptr, "HE") == 0) {
    M = MOLE_WT_HE;
    atoms = GAUSSIAN_MONATOMIC;
    gas = GAS_HE;
    omega = OMEGA_HE;
    mu_not = MU_NOT_HE;
  } else if (strcmp(string_ptr, "N2") == 0) {
    M = MOLE_WT_N2;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_N2;
    omega = OMEGA_N2;
    mu_not = MU_NOT_N2;
  } else if (strcmp(string_ptr, "O2") == 0) {
    M = MOLE_WT_O2;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_O2;
    omega = OMEGA_O2;
    mu_not = MU_NOT_O2;
  } else {
    cout << "*****URECOGNIZED GAS TYPE!  USING AIR!!!***********\n";
    M = MOLE_WT_AIR;
    atoms = GAUSSIAN_DIATOMIC;
    gas = GAS_AIR;
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
  } /* endif */
  if(atoms == GAUSSIAN_MONATOMIC){
    erot = 0.0;
  }
}

/*********************************************************
 * Gaussian2D_pState::pressure -- Thermodynamic pressure.*
 *********************************************************/

inline double Gaussian2D_pState::pressure(void) {
  assert( atoms == GAUSSIAN_MONATOMIC || atoms == GAUSSIAN_DIATOMIC );
  if(atoms==GAUSSIAN_MONATOMIC){
    return (p.xx+p.yy+p.zz)/3.0;
  }else{
    return (p.xx+p.yy+p.zz+2.0*erot)/5.0;
  }
}

inline double Gaussian2D_pState::pressure(void) const {
  assert( atoms == GAUSSIAN_MONATOMIC || atoms == GAUSSIAN_DIATOMIC );
  if(atoms==GAUSSIAN_MONATOMIC){
    return (p.xx+p.yy+p.zz)/3.0;
  }else{
    return (p.xx+p.yy+p.zz+2.0*erot)/5.0;
  }
}

/*********************************************************
 * Gaussian2D_pState::DetP -- Determinant of Press tensor*
 *********************************************************/

inline double Gaussian2D_pState::DetP(void) {
  return (p.xx*p.yy-sqr(p.xy));
}

inline double Gaussian2D_pState::DetP(void) const {
  return (p.xx*p.yy-sqr(p.xy));
}

/************************************************************
 * Gaussian2D_pState::invalid, checks for physical validity *
 ************************************************************/

inline int Gaussian2D_pState::invalid() const{

  int check(0);
  double det = DetP();

  //ensure density > 0
  //and all relevant pressures > 0

  if(   d <= 0.0   ||
        det <= 0.0 ||
	p.xx <= 0.0 ||
	p.yy <= 0.0 ||
        p.zz <= 0.0 ){
    check = 1;
  }

  if((atoms == GAUSSIAN_DIATOMIC) &&
     (erot <= 0.0))
    { check = 1; }

  return check;
}

inline void Gaussian2D_pState::make_valid() {

  if(p.xx<0.0){p.xx = TOLER;}
  if(p.yy<0.0){p.yy = TOLER;}

  if(DetP() < TOLER/MILLION*(p.xx*p.yy)) {
    p.xy = sgn(p.xy)*(1.0-TOLER)*sqrt(p.xx*p.yy);
  }
  return;
}


/********************************************************
 * Gaussian2D_pState::T -- Temperature.                 *
 ********************************************************/

inline double Gaussian2D_pState::T(void) {
  if(pressure()<=0.0) {
    cout << "P=" << pressure() << endl;
    return(THOUSAND);
  }
  if(d<=0.0) {
    cout << "d=" << d << endl;
    return (THOUSAND);
  }
    return (pressure()/(d*AVOGADRO*THOUSAND/M*BOLTZMANN));
}

inline double Gaussian2D_pState::T(void) const {
  if(pressure()<=0.0) {
    cout << "P=" << pressure() << endl;
    return(THOUSAND);
  }
  if(d<=0.0) {
    cout << "d=" << d << endl;
    return (THOUSAND);
  }
    return (pressure()/(d*AVOGADRO*THOUSAND/M*BOLTZMANN));
}

/********************************************************
 * Gaussian2D_pState::E -- Total energy.                   *
 ********************************************************/

inline Tensor2D Gaussian2D_pState::E(void) {

    Tensor2D temp;

    temp.xx = p.xx + d*v.x*v.x;
    temp.xy = p.xy + d*v.x*v.y;
    temp.yy = p.yy + d*v.y*v.y;
    temp.zz = p.zz;

    return (temp); 
}

inline Tensor2D Gaussian2D_pState::E(void) const {

  Tensor2D temp;
  
  temp.xx = p.xx + d*v.x*v.x;
  temp.xy = p.xy + d*v.x*v.y;
  temp.yy = p.yy + d*v.y*v.y;
  temp.zz = p.zz;
  
  return (temp); 

}

/********************************************************
 * Gaussian2D_pState::h -- Specific enthalpy.              *
 ********************************************************/

inline Tensor2D Gaussian2D_pState::h(void) {

    double hxx, hxy, hyy, hzz;

    hxx = 1.5*p.xx/d+0.5*v.x*v.x;
    hxy = 1.5*p.xy/d+0.5*v.x*v.y;
    hyy = 1.5*p.yy/d+0.5*v.y*v.y;
    hzz = 1.5*p.zz/d;

    return (Tensor2D(hxx,hxy,hyy,hzz));
}

inline Tensor2D Gaussian2D_pState::h(void) const {

    double hxx, hxy, hyy, hzz;

    hxx = 1.5*p.xx/d+0.5*v.x*v.x;
    hxy = 1.5*p.xy/d+0.5*v.x*v.y;
    hyy = 1.5*p.yy/d+0.5*v.y*v.y;
    hzz = 1.5*p.zz/d;

    return (Tensor2D(hxx,hxy,hyy,hzz));
}

/********************************************************
 * Gaussian2D_pState::axx & ayy -- Sound speeds. (sort of) *
 ********************************************************/

inline double Gaussian2D_pState::axx(void) {
  assert( p.xx>0.0 && d>0.0 );
  return (sqrt(p.xx/d));
}

inline double Gaussian2D_pState::axx(void) const {
  assert( p.xx>0.0 && d>0.0 );
  return (sqrt(p.xx/d));
}

inline double Gaussian2D_pState::ayy(void) {
    assert( p.yy>0.0 && d>0.0 );
    return (sqrt(p.yy/d));
}

inline double Gaussian2D_pState::ayy(void) const {
  assert( p.yy>0.0 && d>0.0 );
  return (sqrt(p.yy/d));
}

inline double Gaussian2D_pState::sound(void) {
  if(atoms == 1){
    return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3.0/d));
  }else if (atoms == 2){
    return (sqrt(1.4*(p.xx+p.yy+p.zz+2.0*erot)/5.0/d));
  }else{
    cout << "Error....atoms value of " << atoms << "is not allowed.\n";
    cout << "assuming atoms = 1\n";
    return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3/d));
  }
}

inline double Gaussian2D_pState::sound(void) const {
    if(atoms == 1){
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3.0/d));
    }else if (atoms == 2){
      return (sqrt(1.4*(p.xx+p.yy+p.zz+2.0*erot)/5.0/d));
    }else{
      cout << "Error....atoms value of " << atoms << "is not allowed.\n";
      cout << "assuming atoms = 1\n";
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3/d));
    }
}

/********************************************************
 * Gaussian2D_pState::R -- Specific Gas Constant.       *
 ********************************************************/
inline double Gaussian2D_pState::R(void) const {
  return R_UNIVERSAL/M;
}

/********************************************************
 * Gaussian2D_pState::Gamma                             *
 ********************************************************/
inline double Gaussian2D_pState::Gamma(void) const {
  switch(atoms) {
    case 1 :
      return (5.0/3.0);
    case 2 :
      return 1.4;
    default :
      cout << "Error determining Gamma." << endl;
      return -1.0;
  };
}

/********************************************************
 * Gaussian2D_pState::Cv -- specific heat at constant V.*
 ********************************************************/
inline double Gaussian2D_pState::Cv(void) const {
  return R()/(Gamma()-1.0);
}

/********************************************************
 * Gaussian2D_pState::Cp -- specific heat at constant P.*
 ********************************************************/
inline double Gaussian2D_pState::Cp(void) const {
  double g(Gamma());
  return R()*g/(g-1.0);
}

/********************************************************
 * Gaussian2D_pState::K -- Thermal conductivity.        *
 ********************************************************/
inline double Gaussian2D_pState::K(void) const {
  return Cp()*viscosity()/pr;
}



/********************************************************
 * Gaussian2D_pState::dv -- Momentum.                   *
 ********************************************************/
inline Vector2D Gaussian2D_pState::dv(void) {
  return (d*v);
}

inline Vector2D Gaussian2D_pState::dv(void) const {
  return (d*v);
}

inline double Gaussian2D_pState::dv(const Vector2D &n) {
  return (d*(v*n));
}

inline double Gaussian2D_pState::dv(const Vector2D &n) const {
  return (d*(v*n));
}

/********************************************************
 * Gaussian2D_pState::viscosity -- Viscosity.           *
 ********************************************************/
inline double Gaussian2D_pState::viscosity(void){
  
//  double omega, mu_not;  //These are static variables now
//
//  switch(gas) {
//  case GAS_AIR:
//    omega = OMEGA_AIR;
//    mu_not = MU_NOT_AIR;
//    break;
//  case GAS_A:
//    omega = OMEGA_A;
//    mu_not = MU_NOT_A;
//    break;
//  case GAS_CO:
//    omega = OMEGA_CO;
//    mu_not = MU_NOT_CO;
//    break;
//  case GAS_H2:
//    omega = OMEGA_H2;
//    mu_not = MU_NOT_H2;
//    break;
//  case GAS_HE:
//    omega = OMEGA_HE;
//    mu_not = MU_NOT_HE;
//    break;
//  case GAS_N2:
//    omega = OMEGA_N2;
//    mu_not = MU_NOT_N2;
//    break;
//  case GAS_O2:
//    omega = OMEGA_O2;
//    mu_not = MU_NOT_O2;
//    break;
//  default:
//    omega = OMEGA_AIR;
//    mu_not = MU_NOT_AIR;
//    break;
//  }  

  return (mu_not*pow((T()/273.0),omega));

}

inline double Gaussian2D_pState::viscosity(void) const{
  
//  double omega, mu_not;  //These are static variables now
//
//  switch(gas) {
//  case GAS_AIR:
//    omega = OMEGA_AIR;
//    mu_not = MU_NOT_AIR;
//    break;
//  case GAS_A:
//    omega = OMEGA_A;
//    mu_not = MU_NOT_A;
//    break;
//  case GAS_CO:
//    omega = OMEGA_CO;
//    mu_not = MU_NOT_CO;
//    break;
//  case GAS_H2:
//    omega = OMEGA_H2;
//    mu_not = MU_NOT_H2;
//    break;
//  case GAS_HE:
//    omega = OMEGA_HE;
//    mu_not = MU_NOT_HE;
//    break;
//  case GAS_N2:
//    omega = OMEGA_N2;
//    mu_not = MU_NOT_N2;
//    break;
//  case GAS_O2:
//    omega = OMEGA_O2;
//    mu_not = MU_NOT_O2;
//    break;
//  default:
//    omega = OMEGA_AIR;
//    mu_not = MU_NOT_AIR;
//    break;
//  }  

  return (mu_not*pow((T()/273.0),omega));

}

inline double Gaussian2D_pState::nu(void){
  return (viscosity()/d);
}

inline double Gaussian2D_pState::nu(void) const{
  return (viscosity()/d);
}

inline double Gaussian2D_pState::bulk_viscosity(void){
  return (3.0*viscosity());
}

inline double Gaussian2D_pState::bulk_viscosity(void) const{
  return (3.0*viscosity());
}

/********************************************************
 * Gaussian2D_pState -- Mean Free Path (hard spheres)   *
 ********************************************************/
inline double Gaussian2D_pState::mfp(void){
  if(pressure() < TOLER || d < TOLER) {
    return MILLION;
  } else {
    return (16.0*viscosity()/(5.0*sqrt(2.0*PI*d*pressure())));
  }
}

/**********************************************************************
 * Gaussian2D_pState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Gaussian2D_pState::burningrate(void) const {
  return -BETA_APHTPB*pow(pressure(),N_APHTPB);
}

/**********************************************************************
 * Gaussian2D_pState -- Assignment operator.                          *
 **********************************************************************/
inline Gaussian2D_pState& Gaussian2D_pState::operator =(const Gaussian2D_pState &W) {
  d = W.d; v = W.v; p = W.p; erot = W.erot; 
#ifdef _GAUSSIAN_HEAT_TRANSFER_
  q = W.q;
#endif
  return *this;
}

/********************************************************
 * Gaussian2D_pState -- Binary arithmetic operators.    *
 ********************************************************/
inline Gaussian2D_pState operator +(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (Gaussian2D_pState(W1.d+W2.d,W1.v+W2.v,W1.p+W2.p,W1.erot+W2.erot));
}

inline Gaussian2D_pState operator -(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (Gaussian2D_pState(W1.d-W2.d,W1.v-W2.v,W1.p-W2.p,W1.erot-W2.erot));
}

// Inner product operator.
inline double operator *(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.p.xx*W2.p.xx+W1.p.xy*W2.p.xy+W1.p.yy*W2.p.yy+
           W1.p.zz*W2.p.zz+W1.erot*W2.erot);
}

inline Gaussian2D_pState operator *(const Gaussian2D_pState &W, const double &a) {
  return (Gaussian2D_pState(a*W.d,a*W.v,a*W.p,a*W.erot));
}

inline Gaussian2D_pState operator *(const double &a, const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(a*W.d,a*W.v,a*W.p,a*W.erot));
}

inline Gaussian2D_pState operator /(const Gaussian2D_pState &W, const double &a) {
  return (Gaussian2D_pState(W.d/a,W.v/a,W.p/a,W.erot/a));
}

// My useful solution state product operator.
inline Gaussian2D_pState operator ^(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
   return (Gaussian2D_pState(W1.d*W2.d,W1.v.x*W2.v.x,W1.v.y*W2.v.y,W1.p.xx*W2.p.xx,
                          W1.p.xy*W2.p.xy,W1.p.yy*W2.p.yy,W1.p.zz*W2.p.zz,W1.erot*W2.erot));
}

/********************************************************
 * Gaussian2D_pState -- Unary arithmetic operators.     *
 ********************************************************/
inline Gaussian2D_pState operator +(const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(W.d,W.v,W.p,W.erot));
}

inline Gaussian2D_pState operator -(const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(-W.d,-W.v,-W.p,-W.erot));
}

/********************************************************
 * Gaussian2D_pState -- Shortcut arithmetic operators.  *
 ********************************************************/
inline Gaussian2D_pState &operator +=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  W1.d    += W2.d;
  W1.v    += W2.v;
  W1.p    += W2.p;
  W1.erot += W2.erot;
  return (W1);
}

inline Gaussian2D_pState &operator -=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  W1.d    -= W2.d;
  W1.v    -= W2.v;
  W1.p    -= W2.p;
  W1.erot -= W2.erot;
  return (W1);
}

inline Gaussian2D_pState& Gaussian2D_pState::operator *=(const double &a) {
  d *= a; v.x *= a; v.y *= a; p.xx *= a;
  p.xy *=a; p.yy *= a; p.zz *= a; erot *=a;
  return *this;
}

inline Gaussian2D_pState& Gaussian2D_pState::operator /=(const double &a) {
  d /= a; v.x /= a; v.y /= a; p.xx /= a;
  p.xy /=a; p.yy /= a; p.zz /= a; erot /=a;
  return *this;
}


/********************************************************
 * Gaussian2D_pState -- Relational operators.           *
 ********************************************************/
inline int operator ==(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.p == W2.p && W1.erot == W2.erot);
}

inline int operator !=(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.p != W2.p || W1.erot != W2.erot);
}

/********************************************************
 * Gaussian2D_pState -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Gaussian2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.d  << " " << W.v.x << " " << W.v.y << " " << W.p.xx
           << " " << W.p.xy << " " << W.p.yy << " " << W.p.zz << " " << W.erot;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Gaussian2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v.x >> W.v.y >> W.p.xx >> W.p.xy
          >> W.p.yy >> W.p.zz >> W.erot;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Gaussian2D_cState::setgas -- Assign gas constants.   *
 ********************************************************/
inline void Gaussian2D_cState::setgas(void) {
  M = MOLE_WT_AIR;
  atoms = GAUSSIAN_DIATOMIC;
  gas = GAS_AIR;
}

inline void Gaussian2D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
    } else if (strcmp(string_ptr, "A") == 0) {
      M = MOLE_WT_A;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_A;
   } else if (strcmp(string_ptr, "CO") == 0) {
     M = MOLE_WT_CO;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_CO;
   } else if (strcmp(string_ptr, "H2") == 0) {
     M = MOLE_WT_H2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_H2;
   } else if (strcmp(string_ptr, "HE") == 0) {
     M = MOLE_WT_HE;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_HE;
   } else if (strcmp(string_ptr, "N2") == 0) {
     M = MOLE_WT_N2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_N2;
   } else if (strcmp(string_ptr, "O2") == 0) {
     M = MOLE_WT_O2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_O2;
   } else {
     cout << "*****URECOGNIZED GAS TYPE!  USING AIR!!!***********\n";
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
   } /* endif */
   if(atoms == GAUSSIAN_MONATOMIC){
     erot = 0.0;
   }
}

/********************************************************
 * Gaussian2D_cState::v -- Flow velocity.               *
 ********************************************************/
inline Vector2D Gaussian2D_cState::v(void) {
  return (dv/d);
}

inline Vector2D Gaussian2D_cState::v(void) const {
  return (dv/d);
}

inline double Gaussian2D_cState::v(const Vector2D &n) {
  return ((dv*n)/d);
}

inline double Gaussian2D_cState::v(const Vector2D &n) const {
  return ((dv*n)/d);
}

/********************************************************
 * Gaussian2D_cState::p -- Pressure.                    *
 ********************************************************/

inline Tensor2D Gaussian2D_cState::p(void) {

    Tensor2D temp;

    temp.xx = E.xx - dv.x*dv.x/d;
    temp.xy = E.xy - dv.x*dv.y/d;
    temp.yy = E.yy - dv.y*dv.y/d;
    temp.zz = E.zz;

    return (temp); 

}

inline Tensor2D Gaussian2D_cState::p(void) const {

    Tensor2D temp;

    temp.xx = E.xx - dv.x*dv.x/d;
    temp.xy = E.xy - dv.x*dv.y/d;
    temp.yy = E.yy - dv.y*dv.y/d;
    temp.zz = E.zz;

    return (temp); 

}

/********************************************************
 * Gaussian2D_cState::p -- Pressure.                    *
 ********************************************************/
inline double Gaussian2D_cState::pressure(void) const {

  Tensor2D pres = p();

  assert( atoms == GAUSSIAN_MONATOMIC || atoms == GAUSSIAN_DIATOMIC );
  if(atoms==GAUSSIAN_MONATOMIC){
    return (pres.xx+pres.yy+pres.zz)/3.0;
  }else{
    return (pres.xx+pres.yy+pres.zz+2.0*erot)/5.0;
  }
}

/************************************************************
 * Gaussian2D_cState::invalid, checks for physical validity *
 ************************************************************/

inline int Gaussian2D_cState::invalid() const
{
  int check(0);
  double det;
  Tensor2D pressure;

  pressure = p();

  //find determinant of P tensor

  det = (pressure.xx*pressure.yy-sqr(pressure.xy));

  //ensure density > 0
  //and all relevant pressures > 0

  if(   d <= 0.0   ||
        det <= 0.0 ||
	pressure.xx <= 0.0 ||
	pressure.yy <= 0.0 ||
        pressure.zz <= 0.0 ){
    check = 1;
  }

  if((atoms == GAUSSIAN_DIATOMIC) &&
     (erot <= 0.0))
    { check = 1; }

  return check;
}

/**********************************************************************
 * Gaussian2D_cState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Gaussian2D_cState::burningrate(void) const {
  Gaussian2D_pState temp = W(); //UGLY!...but works for now...and I don't use this
  return -BETA_APHTPB*pow(temp.pressure(),N_APHTPB);
}

/**********************************************************************
 * Gaussian2D_cState -- Assignment operator.                          *
 **********************************************************************/
inline Gaussian2D_cState& Gaussian2D_cState::operator =(const Gaussian2D_cState &U) {
  d = U.d; dv = U.dv; E = U.E; erot = U.erot; 
#ifdef _GAUSSIAN_HEAT_TRANSFER_
  q = U.q;
#endif
  return *this;
}

/********************************************************
 * Gaussian2D_cState -- Binary arithmetic operators.       *
 ********************************************************/

inline Gaussian2D_cState operator +(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E,U1.erot+U2.erot));
}

inline Gaussian2D_cState operator -(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E,U1.erot-U2.erot));
}

// Inner product operator.
inline double operator *(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E.xx*U2.E.xx+U1.E.xy*U2.E.xy+U1.E.yy*U2.E.yy
           +U1.E.zz*U2.E.zz+U1.erot*U2.erot);
}

inline Gaussian2D_cState operator *(const Gaussian2D_cState &U, const double &a) {
  return (Gaussian2D_cState(a*U.d,a*U.dv,a*U.E,a*U.erot));
}

inline Gaussian2D_cState operator *(const double &a, const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(a*U.d,a*U.dv,a*U.E,a*U.erot));
}

inline Gaussian2D_cState operator /(const Gaussian2D_cState &U, const double &a) {
  return (Gaussian2D_cState(U.d/a,U.dv/a,U.E/a,U.erot/a));
}

// My useful solution state product operator.
inline Gaussian2D_cState operator ^(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d*U2.d,U1.dv.x*U2.dv.x,U1.dv.y*U2.dv.y,U1.E.xx*U2.E.xx,
                             U1.E.xy*U2.E.xy,U1.E.yy*U2.E.yy,U1.E.zz*U2.E.zz,U1.erot*U2.erot));
}

/********************************************************
 * Gaussian2D_cState -- Unary arithmetic operators.        *
 ********************************************************/
inline Gaussian2D_cState operator +(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.d,U.dv,U.E,U.erot));
}

inline Gaussian2D_cState operator -(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(-U.d,-U.dv,-U.E,-U.erot));
}

/********************************************************
 * Gaussian2D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Gaussian2D_cState &operator +=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  U1.d    += U2.d;
  U1.dv   += U2.dv;
  U1.E    += U2.E;
  U1.erot += U2.erot;
  return (U1);
}

inline Gaussian2D_cState &operator -=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  U1.d    -= U2.d;
  U1.dv   -= U2.dv;
  U1.E    -= U2.E;
  U1.erot -= U2.erot;
  return (U1);
}

inline Gaussian2D_cState& Gaussian2D_cState::operator *=(const double &a) {
  d *= a; dv.x *= a; dv.y *= a; E.xx *= a;
  E.xy *= a; E.yy *= a; E.zz *= a; erot *= a;
  return *this;
}

inline Gaussian2D_cState& Gaussian2D_cState::operator /=(const double &a) {
  d /= a; dv.x /= a; dv.y /= a; E.xx /= a;
  E.xy /= a; E.yy /= a; E.zz /= a; erot /= a;
  return *this;
}

/********************************************************
 * Gaussian2D_cState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E && U1.erot == U2.erot);
}

inline int operator !=(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E || U1.erot != U2.erot);
}

/********************************************************
 * Gaussian2D_cState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Gaussian2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " << U.E.xx
           << " " << U.E.xy << " " << U.E.yy << " " << U.E.zz << " " << U.erot;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Gaussian2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv.x >> U.dv.y >> U.E.xx >> U.E.xy 
          >> U.E.yy >> U.E.zz >> U.erot;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Gaussian2D_pState::Gaussian2D_pState -- Constructor.    *
 ********************************************************/

inline Gaussian2D_pState::Gaussian2D_pState(const Gaussian2D_cState &U) {
  d = U.d; v = U.v(); p = U.p(); erot = U.erot;
}

/********************************************************
 * Gaussian2D_pState::U -- Conserved solution state.    *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::U(void) {
  return (Gaussian2D_cState(d, dv(), E(), erot));
}

inline Gaussian2D_cState Gaussian2D_pState::U(void) const {
  return (Gaussian2D_cState(d, dv(), E(), erot));
}

inline Gaussian2D_cState Gaussian2D_pState::U(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d, W.dv(), W.E(), W.erot));
}

inline Gaussian2D_cState U(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d, W.dv(), W.E(), W.erot));
}

/********************************************************
 * Gaussian2D_pState::F -- Solution flux (x-direction). *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::F(void) {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::F(void) const {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::F(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

inline Gaussian2D_cState F(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

inline Gaussian2D_cState Gaussian2D_pState::F(const Vector2D &V) const {
  return Gaussian2D_cState(d*(v.x-V.x), 
			   d*(v.x-V.x)*v.x + p.xx, 
			   d*(v.x-V.x)*v.y + p.xy, 
			   d*(v.x-V.x)*v.x*v.x + (v.x-V.x)*p.xx + 2.0*v.x*p.xx,
			   d*(v.x-V.x)*v.x*v.y + (v.x-V.x)*p.xy + v.x*p.xy + v.y*p.xx,
			   d*(v.x-V.x)*v.y*v.y + (v.x-V.x)*p.yy + 2.0*v.y*p.xy,
			   (v.x-V.x)*p.zz,
			   (v.x-V.x)*erot);
}

/********************************************************
 * Gaussian2D_pState::Fx -- Solution flux (x-direction).*
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::Fx(void) {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fx(void) const {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fx(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

inline Gaussian2D_cState Fx(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

/********************************************************
 * Gaussian2D_pState::Fy -- Solution flux (y-direction).   *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::Fy(void) {
  return (Gaussian2D_cState(d*v.y, 
                            d*v.x*v.y + p.xy, 
                            d*v.y*v.y + p.yy, 
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + 2.0*v.y*p.xy + v.x*p.yy,
                            d*v.y*v.y*v.y + 3.0*v.y*p.yy,
                            v.y*p.zz,
                            v.y*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fy(void) const {
  return (Gaussian2D_cState(d*v.y, 
                            d*v.x*v.y + p.xy, 
                            d*v.y*v.y + p.yy, 
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + 2.0*v.y*p.xy + v.x*p.yy,
                            d*v.y*v.y*v.y + 3.0*v.y*p.yy,
                            v.y*p.zz,
                            v.y*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fy(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.y, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.y*W.v.y + W.p.yy, 
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + 2.0*W.v.y*W.p.xy + W.v.x*W.p.yy,
                            W.d*W.v.y*W.v.y*W.v.y + 3.0*W.v.y*W.p.yy,
                            W.v.y*W.p.zz,
                            W.v.y*W.erot));
}

inline Gaussian2D_cState Fy(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.y, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.y*W.v.y + W.p.yy, 
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + 2.0*W.v.y*W.p.xy + W.v.x*W.p.yy,
                            W.d*W.v.y*W.v.y*W.v.y + 3.0*W.v.y*W.p.yy,
                            W.v.y*W.p.zz,
                            W.v.y*W.erot));
}

/********************************************************
 * Gaussian2D_pState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
/*
inline Gaussian2D_cState Gaussian2D_pState::Fn(void) {
  return (Gaussian2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Gaussian2D_cState Gaussian2D_pState::Fn(void) const {
  return (Gaussian2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Gaussian2D_cState Gaussian2D_pState::Fn(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Gaussian2D_cState Fn(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}
*/

#ifdef _GAUSSIAN_HEAT_TRANSFER_
/********************************************************
 * Gaussian2D_pState::Gx & Gy -- elliptic (heat) flux.  *
 ********************************************************/
inline Gaussian2D_cState Gaussian2D_pState::Gx(const Gaussian2D_pState &dWdx) const {

  return Gaussian2D_cState(ZERO,ZERO,ZERO,
			   -q.xxx,
			   -q.xxy,
			   -q.xyy,
			   -q.xzz,
			   ZERO);
}

inline Gaussian2D_cState Gaussian2D_pState::Gy(const Gaussian2D_pState &dWdy) const {

  return Gaussian2D_cState(ZERO,ZERO,ZERO,
			   -q.xxy,
			   -q.xyy,
			   -q.yyy,
			   -q.yzz,
			   ZERO);
}

/********************************************************
 * Gaussian2D_pState::ComputeHeatTerms                  *
 ********************************************************/
inline void Gaussian2D_pState::ComputeHeatTerms(const Gaussian2D_pState &dWdx,
						const Gaussian2D_pState &dWdy,
						const Vector2D &X,//this variable is added to be consistend with Jai
						const int &Axisymmetric) {

  double tau = pr*tt();

  q.xxx = -tau*3.0*(p.xx*(dWdx.p.xx-p.xx/d*dWdx.d)/d +
		    p.xy*(dWdy.p.xx-p.xx/d*dWdy.d)/d);

  q.xxy = -tau*(2.0*(p.xx*(dWdx.p.xy-p.xy/d*dWdx.d)/d +
		     p.xy*(dWdy.p.xy-p.xy/d*dWdy.d)/d) +
	        (p.xy*(dWdx.p.xx-p.xx/d*dWdx.d)/d +
		 p.yy*(dWdy.p.xx-p.xx/d*dWdy.d)/d));

  q.xyy = -tau*(2.0*(p.xy*(dWdx.p.xy-p.xy/d*dWdx.d)/d +
		     p.yy*(dWdy.p.xy-p.xy/d*dWdy.d)/d) +
	        (p.xx*(dWdx.p.yy-p.yy/d*dWdx.d)/d +
		 p.xy*(dWdy.p.yy-p.yy/d*dWdy.d)/d));

  q.yyy = -tau*3.0*(p.xy*(dWdx.p.yy-p.yy/d*dWdx.d)/d +
		    p.yy*(dWdy.p.yy-p.yy/d*dWdy.d)/d);

  q.xzz = -tau*(p.xx*(dWdx.p.zz-p.zz/d*dWdx.d)/d +
		p.xy*(dWdy.p.zz-p.zz/d*dWdy.d)/d);

  q.yzz = -tau*(p.xy*(dWdx.p.zz-p.zz/d*dWdx.d)/d +
		p.yy*(dWdy.p.zz-p.zz/d*dWdy.d)/d);

  return;
}
#endif

/************************************************************
 * Gaussian2D_pState::lambda -- Eigenvalue(s) (x-direction).*
 ************************************************************/

inline Gaussian2D_pState Gaussian2D_pState::lambda(void) {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda(void) const {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double Gaussian2D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double lambda(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.x-sqrt(3.0)*W.axx());
    case 2 :
      return (W.v.x-W.axx());
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x);
    case 5 :
      return (W.v.x);
    case 6 :
      return (W.v.x);
    case 7 :
      return (W.v.x+W.axx());
    case 8 :
      return (W.v.x+sqrt(3.0)*W.axx());
    default:
      return (W.v.x);
  };
}

/************************************************************
 * Gaussian2D_pState::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lambda_x(void) {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_x(void) const {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_x(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda_x(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda_x(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double Gaussian2D_pState::lambda_x(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double lambda_x(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.x-sqrt(3.0)*W.axx());
    case 2 :
      return (W.v.x-W.axx());
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x);
    case 5 :
      return (W.v.x);
    case 6 :
      return (W.v.x);
    case 7 :
      return (W.v.x+W.axx());
    case 8 :
      return (W.v.x+sqrt(3.0)*W.axx());
    default:
      return (W.v.x);
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_x(const Vector2D &V) const {
  double c = axx();
  return (Gaussian2D_pState(v.x-V.x-sqrt(3.0)*c, 
			    v.x-V.x-c, 
			    v.x-V.x, 
			    v.x-V.x, 
			    v.x-V.x, 
			    v.x-V.x, 
			    v.x-V.x+c, 
			    v.x-V.x+sqrt(3.0)*c));
}
/************************************************************
 * Gaussian2D_pState::lambda_y -- Eigenvalue(s) (y-direction). *
 ************************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lambda_y(void) {
  double c = ayy();
  return (Gaussian2D_pState(v.y-sqrt(3.0)*c, v.y-c, v.y, v.y, v.y, v.y, v.y+c, v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_y(void) const {
  double c = ayy();
  return (Gaussian2D_pState(v.y-sqrt(3.0)*c, v.y-c, v.y, v.y, v.y, v.y, v.y+c, v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_y(const Gaussian2D_pState &W) {
  double c = W.ayy();
  return (Gaussian2D_pState(W.v.y-sqrt(3.0)*c, W.v.y-c, W.v.y, W.v.y, W.v.y, W.v.y, W.v.y+c, W.v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda_y(const Gaussian2D_pState &W) {
  double c = W.ayy();
  return (Gaussian2D_pState(W.v.y-sqrt(3.0)*c, W.v.y-c, W.v.y, W.v.y, W.v.y, W.v.y, W.v.y+c, W.v.y+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda_y(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.y-sqrt(3.0)*ayy());
    case 2 :
      return (v.y-ayy());
    case 3 :
      return (v.y);
    case 4 :
      return (v.y);
    case 5 :
      return (v.y);
    case 6 :
      return (v.y);
    case 7 :
      return (v.y+ayy());
    case 8 :
      return (v.y+sqrt(3.0)*ayy());
    default:
      return (v.y);
  };
}

inline double Gaussian2D_pState::lambda_y(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.y-sqrt(3.0)*ayy());
    case 2 :
      return (v.y-ayy());
    case 3 :
      return (v.y);
    case 4 :
      return (v.y);
    case 5 :
      return (v.y);
    case 6 :
      return (v.y);
    case 7 :
      return (v.y+ayy());
    case 8 :
      return (v.y+sqrt(3.0)*ayy());
    default:
      return (v.y);
  };
}

inline double lambda_y(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.y-sqrt(3.0)*W.ayy());
    case 2 :
      return (W.v.y-W.ayy());
    case 3 :
      return (W.v.y);
    case 4 :
      return (W.v.y);
    case 5 :
      return (W.v.y);
    case 6 :
      return (W.v.y);
    case 7 :
      return (W.v.y+W.ayy());
    case 8 :
      return (W.v.y+sqrt(3.0)*W.ayy());
    default:
      return (W.v.y);
  };
}

/********************************************************
 * Gaussian2D_pState::rc -- Conserved right eigenvector *
 *                       (x-direction).                 *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::rc(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x-c,2.0*(W.v.y-W.p.xy/(c*W.d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x+c,2.0*(W.v.y+W.p.xy/(c*W.d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*c, W.v.y+sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c+2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c+W.v.x*sqrt(3.0)*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c+2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
  };
}

/*********************************************************
 * Gaussian2D_pState::rc_x -- Conserved right eigenvector*
 *                         (x-direction).                *
 *********************************************************/
inline Gaussian2D_cState Gaussian2D_pState::rc_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc_x(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x-c,2.0*(W.v.y-W.p.xy/(c*W.d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x+c,2.0*(W.v.y+W.p.xy/(c*W.d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*c, W.v.y+sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c+2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c+W.v.x*sqrt(3.0)*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c+2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
  };
}

/*********************************************************
 * Gaussian2D_pState::rc_y -- Conserved right eigenvector*
 *                         (y-direction).                *
 *********************************************************/
inline Gaussian2D_cState Gaussian2D_pState::rc_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x-p.xy/(c*d)),v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x+p.xy/(c*d)),v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*p.xy/(c*d),v.y+sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c+2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c+sqrt(3.0)*v.y*p.xy+sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c+2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x-p.xy/(c*d)),v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x+p.xy/(c*d)),v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*p.xy/(c*d),v.y+sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c+2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c+sqrt(3.0)*v.y*p.xy+sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c+2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc_y(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*W.p.xy/(c*W.d),W.v.y-sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c-2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c-sqrt(3.0)*W.v.y*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c-2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(W.v.x-W.p.xy/(c*W.d)),W.v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(W.v.x+W.p.xy/(c*W.d)),W.v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*W.p.xy/(c*W.d),W.v.y+sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c+2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c+sqrt(3.0)*W.v.y*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c+2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*W.p.xy/(c*W.d),W.v.y-sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c-2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c-sqrt(3.0)*W.v.y*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c-2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
  };
}


/********************************************************
 * Gaussian2D_pState::lp -- Primitive left eigenvector     *
 *                       (x-direction).                 *
 ********************************************************/

inline Gaussian2D_pState Gaussian2D_pState::lp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,W.p.xy/(2.0*W.d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*W.p.xy*W.p.xy-W.d*c*c*W.p.yy)/(3.0*W.d*W.d*c*c*c*c),
                                -2.0*W.p.xy/(c*c*W.d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,-W.p.xy/(2.0*W.d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

/********************************************************
 * Gaussian2D_pState::lp_x -- Primitive left eigenvector   *
 *                         (x-direction).               *
 ********************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp_x(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,W.p.xy/(2.0*W.d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*W.p.xy*W.p.xy-W.d*c*c*W.p.yy)/(3.0*W.d*W.d*c*c*c*c),
                                -2.0*W.p.xy/(c*c*W.d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,-W.p.xy/(2.0*W.d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

/********************************************************
 * Gaussian2D_pState::lp_y -- Primitive left eigenvector*
 *                         (y-direction).               *
 ********************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,-ONE/(2*c),p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*p.xy/(d*c*c),(4.0*p.xy*p.xy-p.xx*d*c*c)/(3.0*d*d*c*c*c*c),
                                ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,ONE/(2*c),-p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,-ONE/(2*c),p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*p.xy/(d*c*c),(4.0*p.xy*p.xy-p.xx*d*c*c)/(3.0*d*d*c*c*c*c),
                                ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,ONE/(2*c),-p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp_y(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,W.d/2.0,-W.p.xy/(2.0*c*c),ZERO,-ONE/(2*c),W.p.xy/(2*W.d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*W.p.xy/(W.d*c*c),
                                (4.0*W.p.xy*W.p.xy-W.p.xx*W.d*c*c)/(3.0*W.d*W.d*c*c*c*c),ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,W.d/2.0,-W.p.xy/(2.0*c*c),ZERO,ONE/(2*c),-W.p.xy/(2*W.d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

/*************************************************************
 * Gaussian2D_pState::S -- Source terms (axisymmetric flow). *
 *************************************************************/
inline Gaussian2D_cState Gaussian2D_pState::S(const Vector2D &X) {
  return (Gaussian2D_cState(-d*v.y/X.y,
			    -(d*v.x*v.y+p.xy)/X.y,
			    -(d*v.y*v.y+p.yy-p.zz)/X.y,
			    -(d*v.x*v.x*v.y+2.0*v.x*p.xy+v.y*p.xx)/X.y,
			    -(-(v.x*p.zz)+(d*v.x*v.y*v.y+2.0*v.y*p.xy+v.x*p.yy))/X.y,
			    -(-2.0*(v.y*p.zz)+(d*v.y*v.y*v.y+3.0*v.y*p.yy))/X.y,
			    -(3.0*(v.y*p.zz))/X.y,  //3?
			    -erot*v.y/X.y));
}

inline Gaussian2D_cState Gaussian2D_pState::S(const Vector2D &X) const {
  return (Gaussian2D_cState(-d*v.y/X.y,
			    -(d*v.x*v.y+p.xy)/X.y,
			    -(d*v.y*v.y+p.yy-p.zz)/X.y,
			    -(d*v.x*v.x*v.y+2.0*v.x*p.xy+v.y*p.xx)/X.y,
			    -(-(v.x*p.zz)+(d*v.x*v.y*v.y+2.0*v.y*p.xy+v.x*p.yy))/X.y,
			    -(-2.0*(v.y*p.zz)+(d*v.y*v.y*v.y+3.0*v.y*p.yy))/X.y,
			    -(3.0*(v.y*p.zz))/X.y,  //3?
			    -erot*v.y/X.y));
}

inline Gaussian2D_cState S(const Gaussian2D_pState &W, const Vector2D &X) {
  return (Gaussian2D_cState(-W.d*W.v.y/X.y,
			    -(W.d*W.v.x*W.v.y+W.p.xy)/X.y,
			    -(W.d*W.v.y*W.v.y+W.p.yy-W.p.zz)/X.y,
			    -(W.d*W.v.x*W.v.x*W.v.y+2.0*W.v.x*W.p.xy+W.v.y*W.p.xx)/X.y,
			    -(-(W.v.x*W.p.zz)+(W.d*W.v.x*W.v.y*W.v.y+2.0*W.v.y*W.p.xy+W.v.x*W.p.yy))/X.y,
			    -(-2.0*(W.v.y*W.p.zz)+(W.d*W.v.y*W.v.y*W.v.y+3.0*W.v.y*W.p.yy))/X.y,
			    -(3.0*(W.v.y*W.p.zz))/X.y,  //3?
			    -W.erot*W.v.y/X.y));
}

/*************************************************************
 * Gaussian2D_pState -- relatiation times                    *
 *************************************************************/
inline double Gaussian2D_pState::tt() const {
  return viscosity()/pressure();
}

inline double Gaussian2D_pState::tr() const {
  return 15.0/4.0*bulk_viscosity()/pressure();
}

/*************************************************************
 * Gaussian2D_pState -- temperature jump distance            *
 *                                                           *
 *      This is used for temperature slip BC.                *
 *                                                           *
 * For more details, see "Kinetic Theory of Gases" by        *
 * Kennard(1938) p.311-314.                                  *
 *      Available at the UTIAS library.                      *
 *                                                           *
 * Kennard simply uses "g", but I don't want it confused     *
 * with gamma.                                               *
 *************************************************************/
inline double Gaussian2D_pState::gt() const {
  return (2.0-alpha_t)/alpha_t * sqrt(2*PI*R()*T()) * K()/((Gamma()+1.0)*Cv()*pressure());
}

/*************************************************************
 * Gaussian2D_pState -- set_temperatur_d                     *
 *************************************************************/
inline void Gaussian2D_pState::set_temperature_d(double temperature) {
  //This uses the "thermodynamic" pressure and does not change
  //velocities, therefore it implicitly changes Mach numbers.
  d = pressure()*M/(AVOGADRO*BOLTZMANN*temperature)/THOUSAND;
  return;
}

/********************************************************
 * Gaussian2D_cState::Gaussian2D_cState -- Constructor. *
 ********************************************************/

inline Gaussian2D_cState::Gaussian2D_cState(const Gaussian2D_pState &W) {
  d = W.d; dv = W.dv(); E = W.E(); erot = W.erot;
}

/********************************************************
 * Gaussian2D_cState::W -- Primitive solution state.    *
 ********************************************************/

inline Gaussian2D_pState Gaussian2D_cState::W(void) {
  return (Gaussian2D_pState(d, v(), p(), erot));
}

inline Gaussian2D_pState Gaussian2D_cState::W(void) const {
  return (Gaussian2D_pState(d, v(), p(), erot));
}

inline Gaussian2D_pState Gaussian2D_cState::W(const Gaussian2D_cState &U) {
  return (Gaussian2D_pState(U.d, U.v(), U.p(), U.erot));
}

inline Gaussian2D_pState W(const Gaussian2D_cState &U) {
  return (Gaussian2D_pState(U.d, U.v(), U.p(), U.erot));
}

/********************************************************
 * Gaussian2D_cState::F -- Solution flux (x-direction). *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::F(void) {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + pressure.xx,
                            dv.x*dv.y/d + pressure.xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*pressure.xx,
                            E.xy*dv.x/d + dv.x/d*pressure.xy + dv.y/d*pressure.xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*pressure.xy,
                            E.zz*dv.x/d,
			    erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::F(void) const {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + pressure.xx,
                            dv.x*dv.y/d + pressure.xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*pressure.xx,
                            E.xy*dv.x/d + dv.x/d*pressure.xy + dv.y/d*pressure.xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*pressure.xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::F(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + pressure.xx,
                            U.dv.x*U.dv.y/U.d + pressure.xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*pressure.xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*pressure.xy + U.dv.y/U.d*pressure.xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*pressure.xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

inline Gaussian2D_cState F(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.x, 
			    U.dv.x*U.dv.x/U.d + pressure.xx,
                            U.dv.x*U.dv.y/U.d + pressure.xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*pressure.xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*pressure.xy + U.dv.y/U.d*pressure.xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*pressure.xy,
                            U.E.zz*U.dv.x/U.d,
			    U.erot*U.dv.x/U.d));
}

inline Gaussian2D_cState Gaussian2D_cState::F(const Vector2D &V) const {
  Tensor2D pressure = p();
  Vector2D vel = v();
  return Gaussian2D_cState(d*(vel.x-V.x), 
			   d*(vel.x-V.x)*vel.x + pressure.xx, 
			   d*(vel.x-V.x)*vel.y + pressure.xy, 
			   d*(vel.x-V.x)*vel.x*vel.x + (vel.x-V.x)*pressure.xx + 2.0*vel.x*pressure.xx,
			   d*(vel.x-V.x)*vel.x*vel.y + (vel.x-V.x)*pressure.xy + vel.x*pressure.xy + vel.y*pressure.xx,
			   d*(vel.x-V.x)*vel.y*vel.y + (vel.x-V.x)*pressure.yy + 2.0*vel.y*pressure.xy,
			   (vel.x-V.x)*pressure.zz,
			   (vel.x-V.x)*erot);
}

/********************************************************
 * Gaussian2D_cState::Fx -- Solution flux (x-direction).   *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::Fx(void) {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + pressure.xx,
                            dv.x*dv.y/d + pressure.xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*pressure.xx,
                            E.xy*dv.x/d + dv.x/d*pressure.xy + dv.y/d*pressure.xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*pressure.xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fx(void) const {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + pressure.xx,
                            dv.x*dv.y/d + pressure.xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*pressure.xx,
                            E.xy*dv.x/d + dv.x/d*pressure.xy + dv.y/d*pressure.xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*pressure.xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fx(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + pressure.xx,
                            U.dv.x*U.dv.y/U.d + pressure.xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*pressure.xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*pressure.xy + U.dv.y/U.d*pressure.xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*pressure.xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

inline Gaussian2D_cState Fx(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + pressure.xx,
                            U.dv.x*U.dv.y/U.d + pressure.xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*pressure.xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*pressure.xy + U.dv.y/U.d*pressure.xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*pressure.xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

/********************************************************
 * Gaussian2D_cState::Fy -- Solution flux (y-direction).*
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::Fy(void) {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.y,
                            dv.y*dv.x/d + pressure.xy,
                            dv.y*dv.y/d + pressure.yy, 
                            E.xx*dv.y/d + 2.0*dv.x/d*pressure.xy,
                            E.xy*dv.y/d + dv.y/d*pressure.xy + dv.x/d*pressure.yy,
                            E.yy*dv.y/d + 2.0*dv.y/d*pressure.yy,
                            E.zz*dv.y/d,
                            erot*dv.y/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fy(void) const {
  Tensor2D pressure = p();
  return (Gaussian2D_cState(dv.y,
                            dv.y*dv.x/d + pressure.xy,
                            dv.y*dv.y/d + pressure.yy, 
                            E.xx*dv.y/d + 2.0*dv.x/d*pressure.xy,
                            E.xy*dv.y/d + dv.y/d*pressure.xy + dv.x/d*pressure.yy,
                            E.yy*dv.y/d + 2.0*dv.y/d*pressure.yy,
                            E.zz*dv.y/d,
                            erot*dv.y/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fy(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.y,
                            U.dv.y*U.dv.x/U.d + pressure.xy,
                            U.dv.y*U.dv.y/U.d + pressure.yy, 
                            U.E.xx*U.dv.y/U.d + 2.0*U.dv.x/U.d*pressure.xy,
                            U.E.xy*U.dv.y/U.d + U.dv.y/U.d*pressure.xy + U.dv.x/U.d*pressure.yy,
                            U.E.yy*U.dv.y/U.d + 2.0*U.dv.y/U.d*pressure.yy,
                            U.E.zz*U.dv.y/U.d,
                            U.erot*U.dv.y/U.d));
}

inline Gaussian2D_cState Fy(const Gaussian2D_cState &U) {
  Tensor2D pressure = U.p();
  return (Gaussian2D_cState(U.dv.y,
                            U.dv.y*U.dv.x/U.d + pressure.xy,
                            U.dv.y*U.dv.y/U.d + pressure.yy, 
                            U.E.xx*U.dv.y/U.d + 2.0*U.dv.x/U.d*pressure.xy,
                            U.E.xy*U.dv.y/U.d + U.dv.y/U.d*pressure.xy + U.dv.x/U.d*pressure.yy,
                            U.E.yy*U.dv.y/U.d + 2.0*U.dv.y/U.d*pressure.yy,
                            U.E.zz*U.dv.y/U.d,
                            U.erot*U.dv.y/U.d));
}

/********************************************************
 * Gaussian2D_cState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
/*
inline Gaussian2D_cState Gaussian2D_cState::Fn(void) {
  return (Gaussian2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fn(void) const {
  return (Gaussian2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fn(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Gaussian2D_cState Fn(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline void Gaussian2D_cState::dFndU(DenseMatrix &dFndU) {
  W().dFndU(dFndU);
}

inline void Gaussian2D_cState::dFndU(DenseMatrix &dFndU) const {
  W().dFndU(dFndU);
}

inline void dFndU(DenseMatrix &dFndU, const Gaussian2D_cState &U) {
  U.W().dFndU(dFndU);
}
*/

#ifdef _GAUSSIAN_HEAT_TRANSFER_
/********************************************************
 * Gaussian2D_cState::Gx & Gy -- elliptic (heat) flux.  *
 ********************************************************/
inline Gaussian2D_cState Gaussian2D_cState::Gx(const Gaussian2D_pState &dWdx) const {

  return Gaussian2D_cState(ZERO,ZERO,ZERO,
			   -q.xxx,
			   -q.xxy,
			   -q.xyy,
			   -q.xzz,
			   ZERO);
}

inline Gaussian2D_cState Gaussian2D_cState::Gy(const Gaussian2D_pState &dWdy) const {

  return Gaussian2D_cState(ZERO,ZERO,ZERO,
			   -q.xxy,
			   -q.xyy,
			   -q.yyy,
			   -q.yzz,
			   ZERO);
}

/********************************************************
 * Gaussian2D_cState::ComputeHeatTerms                  *
 ********************************************************/
inline void Gaussian2D_cState::ComputeHeatTerms(const Gaussian2D_pState &dWdx,
						const Gaussian2D_pState &dWdy,
						const Vector2D &X,//this variable is added to be consistend with Jai
						const int &Axisymmetric) {  //this is very poorly coded.
                                                                            //fix this later
  Tensor2D pressure = p();
  double tau = pr*W().tt();

  q.xxx = -tau*3.0*(pressure.xx*(dWdx.p.xx-pressure.xx/d*dWdx.d)/d +
		    pressure.xy*(dWdy.p.xx-pressure.xx/d*dWdy.d)/d);

  q.xxy = -tau*(2.0*(pressure.xx*(dWdx.p.xy-pressure.xy/d*dWdx.d)/d +
		     pressure.xy*(dWdy.p.xy-pressure.xy/d*dWdy.d)/d) +
	        (pressure.xy*(dWdx.p.xx-pressure.xx/d*dWdx.d)/d +
		 pressure.yy*(dWdy.p.xx-pressure.xx/d*dWdy.d)/d));

  q.xyy = -tau*(2.0*(pressure.xy*(dWdx.p.xy-pressure.xy/d*dWdx.d)/d +
		     pressure.yy*(dWdy.p.xy-pressure.xy/d*dWdy.d)/d) +
	        (pressure.xx*(dWdx.p.yy-pressure.yy/d*dWdx.d)/d +
		 pressure.xy*(dWdy.p.yy-pressure.yy/d*dWdy.d)/d));

  q.yyy = -tau*3.0*(pressure.xy*(dWdx.p.yy-pressure.yy/d*dWdx.d)/d +
		    pressure.yy*(dWdy.p.yy-pressure.yy/d*dWdy.d)/d);

  q.xzz = -tau*(pressure.xx*(dWdx.p.zz-pressure.zz/d*dWdx.d)/d +
		pressure.xy*(dWdy.p.zz-pressure.zz/d*dWdy.d)/d);

  q.yzz = -tau*(pressure.xy*(dWdx.p.zz-pressure.zz/d*dWdx.d)/d +
		pressure.yy*(dWdy.p.zz-pressure.zz/d*dWdy.d)/d);

  return;
}
#endif

/**********************************************************
 * Gaussian2D_cState::S -- Source terms (axisymmetric flow). *
 **********************************************************/
/*
inline Gaussian2D_cState Gaussian2D_cState::S(const Vector2D &X) {
  return (Gaussian2D_cState(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0));
}

inline Gaussian2D_cState Gaussian2D_cState::S(const Vector2D &X) const {
  return (Gaussian2D_cState(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0));
}

inline Gaussian2D_cState S(const Gaussian2D_cState &U, const Vector2D &X) {
  return (Gaussian2D_cState(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0));
}
*/

/********************************************************
 * Useful 2D Gaussian state constants.                  *
 ********************************************************/

const Gaussian2D_pState Gaussian2D_W_STDATM(DENSITY_STDATM,
				      Vector2D_ZERO, PRESSURE_STDATM);
const Gaussian2D_pState Gaussian2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO);
const Gaussian2D_pState Gaussian2D_W_ONE(ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE);
const Gaussian2D_cState Gaussian2D_U_STDATM(Gaussian2D_W_STDATM);
const Gaussian2D_cState Gaussian2D_U_VACUUM(Gaussian2D_W_VACUUM);
const Gaussian2D_cState Gaussian2D_U_ONE(ONE, ONE, ONE, ONE, ONE, ONE, ONE , ONE);

/********************************************************
 * Gaussian2DState -- External subroutines.             *
 ********************************************************/

extern Gaussian2D_pState RoeAverage(const Gaussian2D_pState &Wl,
				    const Gaussian2D_pState &Wr);

extern Gaussian2D_pState Translate(const Gaussian2D_pState &W, const Vector2D &V);

extern Gaussian2D_pState Rotate(const Gaussian2D_pState &W,
	      	                const Vector2D &norm_dir);

extern Gaussian2D_cState Rotate(const Gaussian2D_cState &U,
	      	                const Vector2D &norm_dir);

extern Gaussian2D_pState Reflect(const Gaussian2D_pState &W,
	      	                 const Vector2D &norm_dir);

extern Gaussian2D_pState Reflect(const Gaussian2D_pState &W,
				 const Vector2D &norm_dir,
				 const Vector2D &V);

extern Gaussian2D_pState NoSlip(const Gaussian2D_pState &W,
	      	                const Vector2D &norm_dir);

extern double dTdn(const Gaussian2D_pState &W,
		   const Gaussian2D_pState &dWdx,
		   const Gaussian2D_pState &dWdy,
		   const Vector2D &norm_dir);

extern double Slip_T(const Gaussian2D_pState &W,
		     const double &T,
		     const double &old_T,
		     const Gaussian2D_pState &dWdx,
		     const Gaussian2D_pState &dWdy,
		     const Vector2D &norm_dir);

extern Gaussian2D_pState Adiabatic_Wall(const Gaussian2D_pState &W,
					const Gaussian2D_pState &Wo,
					const Vector2D &norm_dir,
					const Vector2D &V);

extern Gaussian2D_pState Adiabatic_Wall(const Gaussian2D_pState &W,
					const Vector2D &V,
					const Vector2D &norm_dir);

extern Gaussian2D_pState Isothermal_Wall(const Gaussian2D_pState &W,
					 const Vector2D &V,
					 const double &T,
					 const Vector2D &norm_dir);

extern Gaussian2D_pState Isothermal_Wall_Slip_T(const Gaussian2D_pState &W,
						const Vector2D &V,
						const double &T,
						const double &old_T,
						const Gaussian2D_pState &dWdx,
						const Gaussian2D_pState &dWdy,
						const Vector2D &norm_dir);

extern Gaussian2D_pState Knudsen_Layer_Adiabatic(const Gaussian2D_pState &W,
						 const Vector2D &v,
						 const Vector2D &norm_dir);

extern Gaussian2D_pState Knudsen_Layer_Isothermal(const Gaussian2D_pState &W,
						  const Vector2D &V,
						  const double &T,
						  const Vector2D &norm_dir);

extern Gaussian2D_pState Knudsen_Layer_Isothermal_Slip_T(const Gaussian2D_pState &W,
							 const Vector2D &V,
							 const double &T,
							 const double &old_T,
							 const Gaussian2D_pState &dWdx,
							 const Gaussian2D_pState &dWdy,
							 const Vector2D &norm_dir);

extern Gaussian2D_pState RinglebFlowAverageState(const Gaussian2D_pState &Wdum,
						 const Vector2D &Y1,
						 const Vector2D &Y2,
						 const Vector2D &Y3,
						 const Vector2D &Y4);

extern Gaussian2D_pState BC_Characteristic_Pressure(const Gaussian2D_pState &Wi,
                                                    const Gaussian2D_pState &Wo,
	      	                                    const Vector2D &norm_dir);

extern Gaussian2D_pState BC_Characteristic_Velocity(const Gaussian2D_pState &Wi,
                                                    const Gaussian2D_pState &Wo,
						    const Vector2D &norm_dir);

extern Gaussian2D_pState BC_Couette(const Gaussian2D_pState &Wi,
				    const Gaussian2D_pState &Wo);

extern Gaussian2D_pState BC_Developed_Channel_Flow(const Gaussian2D_pState &Wo,
						   const Vector2D &norm_dir,
						   const double &y);

extern Gaussian2D_pState BCs(const Gaussian2D_pState &Wb,
			     const Gaussian2D_pState &Wi,
			     const Gaussian2D_pState &dWdx,
			     const double &dx,
			     const int BC_type,
			     const int End_type);

extern Gaussian2D_pState BCs_x(const Gaussian2D_pState &Wb,
			       const Gaussian2D_pState &Wi,
			       const Gaussian2D_pState &dWdx,
			       const double &dx,
			       const int BC_type,
			       const int End_type);

extern Gaussian2D_pState BCs_y(const Gaussian2D_pState &Wb,
			       const Gaussian2D_pState &Wi,
			       const Gaussian2D_pState &dWdy,
			       const double &dy,
			       const int BC_type,
			       const int End_type);

extern Gaussian2D_pState WaveSpeedPos(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState WaveSpeedNeg(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState WaveSpeedAbs(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixPos(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixNeg(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixAbs(const Gaussian2D_pState &lambda_a,
				      const Gaussian2D_pState &lambda_l,
				      const Gaussian2D_pState &lambda_r);

extern Gaussian2D_cState FluxRoe(const Gaussian2D_pState &Wl,
				 const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe(const Gaussian2D_cState &Ul,
				 const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_x(const Gaussian2D_pState &Wl,
				   const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe_x(const Gaussian2D_cState &Ul,
				   const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_y(const Gaussian2D_pState &Wl,
				   const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe_y(const Gaussian2D_cState &Ul,
				   const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_n(const Gaussian2D_pState &Wl,
				   const Gaussian2D_pState &Wr,
				   const Vector2D &norm_dir);

extern Gaussian2D_cState FluxRoe_n(const Gaussian2D_cState &Ul,
				   const Gaussian2D_cState &Ur,
				   const Vector2D &norm_dir);

extern Gaussian2D_cState FluxRoe_MB(const Gaussian2D_pState &Wl,
				    const Gaussian2D_pState &Wr,
				    const Vector2D &V);

extern Gaussian2D_cState FluxRoe_MB(const Gaussian2D_cState &Ul,
				    const Gaussian2D_cState &Ur,
				    const Vector2D &V);

extern Gaussian2D_cState FluxRoe_MB_n(const Gaussian2D_pState &Wl,
				      const Gaussian2D_pState &Wr,
				      const Vector2D &V,
				      const Vector2D &norm_dir);

extern Gaussian2D_cState FluxRoe_MB_n(const Gaussian2D_cState &Ul,
				      const Gaussian2D_cState &Ur,
				      const Vector2D &V,
				      const Vector2D &norm_dir);

extern Gaussian2D_cState FluxHLLE(const Gaussian2D_pState &Wl,
				  const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE(const Gaussian2D_cState &Ul,
				  const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLE_x(const Gaussian2D_pState &Wl,
				    const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE_x(const Gaussian2D_cState &Ul,
				    const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLE_y(const Gaussian2D_pState &Wl,
				    const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE_y(const Gaussian2D_cState &Ul,
				    const Gaussian2D_cState &Ur);
  
extern Gaussian2D_cState FluxHLLE_n(const Gaussian2D_pState &Wl,
				    const Gaussian2D_pState &Wr,
				    const Vector2D &norm_dir);

extern Gaussian2D_cState FluxHLLE_n(const Gaussian2D_cState &Ul,
				    const Gaussian2D_cState &Ur,
				    const Vector2D &norm_dir);

extern Vector2D HLLE_wavespeeds(const Gaussian2D_pState &Wl,
                                const Gaussian2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Gaussian2D_cState Imposed_adiabatic_wall_n(const Gaussian2D_pState &Wr,
						  const Gaussian2D_pState &Wo,
						  const Vector2D &norm_dir);


/*
  extern Gaussian2D_cState FluxHLLC(const Gaussian2D_pState &Wl,
  const Gaussian2D_pState &Wr);

  extern Gaussian2D_cState FluxHLLC(const Gaussian2D_cState &Ul,
  const Gaussian2D_cState &Ur);

  extern Gaussian2D_cState FluxHLLC_x(const Gaussian2D_pState &Wl,
  const Gaussian2D_pState &Wr);

  extern Gaussian2D_cState FluxHLLC_x(const Gaussian2D_cState &Ul,
  const Gaussian2D_cState &Ur);

  extern Gaussian2D_cState FluxHLLC_y(const Gaussian2D_pState &Wl,
  const Gaussian2D_pState &Wr);

  extern Gaussian2D_cState FluxHLLC_y(const Gaussian2D_cState &Ul,
  const Gaussian2D_cState &Ur);

  extern Gaussian2D_cState FluxHLLC_n(const Gaussian2D_pState &Wl,
  const Gaussian2D_pState &Wr,
  const Vector2D &norm_dir);

  extern Gaussian2D_cState FluxHLLC_n(const Gaussian2D_cState &Ul,
  const Gaussian2D_cState &Ur,
  const Vector2D &norm_dir);
*/

extern Gaussian2D_cState FluxKinetic_x(const Gaussian2D_pState &W1,
				       const Gaussian2D_pState &W2);

extern Gaussian2D_cState HeatFlux_n(const Vector2D &X,
				    Gaussian2D_pState &W,
				    const Gaussian2D_pState &dWdx,
				    const Gaussian2D_pState &dWdy,
				    const Vector2D &norm_dir,
				    const int &Axisymmetric);

extern Gaussian2D_cState HeatFluxDiamondPath_n(const Vector2D &X,
					       const Vector2D &X1, const Gaussian2D_pState &W1,
					       const Vector2D &X2, const Gaussian2D_pState &W2,
					       const Vector2D &X3, const Gaussian2D_pState &W3,
					       const Vector2D &X4, const Gaussian2D_pState &W4,
					       const Vector2D &norm_dir,
					       const int &Axisymmetric,
					       const int &stencil_flag);

extern Gaussian2D_cState HeatFluxHybrid_n(const Vector2D &X,
					  Gaussian2D_pState &W,
					  const Vector2D &X1,
					  const Gaussian2D_pState &W1,
					  const Gaussian2D_pState &dW1dx,
					  const Gaussian2D_pState &dW1dy,
					  const Vector2D &X2,
					  const Gaussian2D_pState &W2,
					  const Gaussian2D_pState &dW2dx,
					  const Gaussian2D_pState &dW2dy,
					  const Vector2D &norm_dir,
					  const int &Axisymmetric);

extern Gaussian2D_pState FlatPlate(const Gaussian2D_pState &Winf,
				   const Vector2D X,
				   double &eta,
				   double &f,
				   double &fp,
				   double &fpp);

extern Gaussian2D_pState Free_Molecular_Exact(const Vector2D position,
					      const Gaussian2D_pState &Winf,
					      Vector2D *nodes,
					      const int number_of_nodes);


extern void Integrate_distribution(Gaussian2D_pState &solution,
				   const double min_angle,
				   const double max_angle,
				   const Gaussian2D_pState &W);



#endif /* _GAUSSIAN2D_STATE_INCLUDED  */

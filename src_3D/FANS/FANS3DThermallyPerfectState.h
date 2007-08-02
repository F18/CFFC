/****************** FANS3DThermallyPerfectState.h *********************
  This class defines the state variables and constructors for the 
  FANS3D thermally perfect gaseous mixture class.
 **********************************************************************/

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED 
#define _FANS3D_THERMALLYPERFECT_STATE_INCLUDED

class FANS3D_ThermallyPerfect_KOmega_cState;
class FANS3D_ThermallyPerfect_KOmega_pState;

// Required C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

// Include required CFFC header files

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../NavierStokes/NavierStokes3DThermallyPerfectState.h"
#endif  //NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

#ifndef  _TURBULENCE_MODELLING_INCLUDED 
#include "../Turbulence/TurbulenceModelling.h"
#endif  //TURBULENCE_MODELLING_INCLUDED

#define NUM_FANS3D_VAR_EXTRA 2  //rho, v(3), p

/********************************************************
 * Class: FANS3D_ThermallyPerfect_KOmega_pState         *
 * (all for thermally perfect mixture                   *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1                           *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      T       -- Return temperature.                  *
 *      e       -- Return total  energy.                *
 *      E       -- Return interal energy.               *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      dv      -- Return momentum.                     *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          *
 *      U       -- Return conserved solution state.     *
 *      Fx      -- Return x-direction solution flux.    *
 *      Fy      -- Return y-direction solution flux.    *
 *      Fz      -- Return z-direction solution flux.    *
 *      Fn      -- Return n-direction solution flux.    *
 *      lambda  -- Return eigenvalue.                   *
 *      rp      -- Return primitive right eigenvector.  *
 *      rc      -- Return conserved right eigenvector.  *
 *      lp      -- Return primitive left eigenvector.   *
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

/*********************************************************


  Density:   rho  kg/m^3
  Velocity:  v    m/s
  Pressure:  p    Pa (N/m^2)
  
  Molecular Mass:          M    kg/mol
  Species Gas Constant:    Rs   J/(kg*K)
  
  Temperature(Kelvin) Dependent data: f(T) 
   
  Heat Capacity (const Pressure):  Cp  J/(kg*K)
  Heat Capacity (const Volume):    Cv  J/(kg*K)
  Specific Heat Ratio:             g
  Specific Enthalpy:               h   J/kg
  Specific Internal Energy:        e   J/kg
  Total Enthalpy:                  H   J/kg 
  Total Internal Energy:           E   J/kg

  Viscosity:                       mu  kg/(m*s) N*s/m^2  
  Thermal Conductivity:            k   N/(s*K)  W.(m*K)

  ns;                number of species
  NASARP1311data *specdata;   Global Species Data

  NUM_VAR_3D;         number of total variables (5+ns)
  Reaction_set React;         Global Reaction Data

  low_temp_range;      low temp data range
  high_temp_range;     high temp data range

***********************************************************/
class FANS3D_ThermallyPerfect_KOmega_pState :public NavierStokes3D_ThermallyPerfect_pState {
  public:
   double k;
   double omega;
     
   Turbulence_Model_k_omega k_omega_model;
 
   
   // constructor
   FANS3D_ThermallyPerfect_KOmega_pState(): NavierStokes3D_ThermallyPerfect_pState(){ 
      k= ONE; omega = MILLION;}
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &value):
      NavierStokes3D_ThermallyPerfect_pState(value){
      k = value; omega = value;}
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const Vector3D &V,
                                         const double &pre, const double &ke, 
                                         const double &Omega) :
      NavierStokes3D_ThermallyPerfect_pState(d, V, pre){
      k = ke; omega = Omega;}

   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const Vector3D &V,
                                         const double &pre, const double &ke, 
                                         const double &Omega, 
                                         const double &frac):
      NavierStokes3D_ThermallyPerfect_pState(d, V, pre, frac){ 
      k = ke; omega = Omega; }
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const double &vx,
                                  const double &vy, const double &vz,
                                  const double &pre, const double &ke, 
                                  const double &Omega,
                                  const double &frac):
      NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac){
      k = ke; omega = Omega;}
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const double &vx,
                                  const double &vy, const double &vz,
                                  const double &ke, const double &Omega,
                                  const double &pre):
      NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre) {
      k = ke; omega = Omega;}
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const double &vx,
                                  const double &vy, const double &vz,
                                  const double &pre, const double &ke, 
                                  const double &Omega, Species *mfrac):
      NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac){
      k = ke; omega = Omega; }
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const Vector3D &V,
                                  const double &pre, const double &ke, 
                                  const double &Omega, Species *mfrac):
      NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){
      k = ke; omega = Omega;  }
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const Vector3D &V,
                                         const double &pre, double *mfrac):
      NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){
      k = ONE; omega = MILLION; 
      set_initial_values(mfrac); }
   
   
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, const Vector3D &V,
                                  const double &pre, const double &ke, 
                                  const double &Omega, double *mfrac):
      NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){
      k = ke; omega = Omega; 
      set_initial_values(mfrac); }
  
      //this is needed for the operator overload returns!!!!
      FANS3D_ThermallyPerfect_KOmega_pState(const  FANS3D_ThermallyPerfect_KOmega_pState &W)
                   {spec = NULL; rho = DENSITY_STDATM; set_initial_values(); Copy(W);}
                    

      static void set_species_data(FANS3D_ThermallyPerfect_KOmega_pState &Wo, 
                                   const int &, 
                                   const string *,
                                   const char *, 
                                   const int &, 
                                   const double&, 
                                   const double *);

   void Vacuum(){
      Euler3D_ThermallyPerfect_pState::Vacuum();
      k = ZERO;
      omega = ZERO;
      
   }

   void Copy(const FANS3D_ThermallyPerfect_KOmega_pState &W);
   double E(void) const ;      //mixture total internal energy
   double H(void) const ;      //mixture total enthalpy
   double Hs(void) const;
   double a(void);                
   double a(void) const;
   
   Tensor3D rotation_tensor(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdy, 
                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   Tensor3D strain_rate(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy, 
                        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   double eddy_viscosity(void) const;      
   double Pr_t(void) const;      
   double Sc_t(void) const;      
   double kappa_t(void) const;      
   double Dm_t(void) const;
  
   // Conserved solution state. /
   FANS3D_ThermallyPerfect_KOmega_cState U(void);
   FANS3D_ThermallyPerfect_KOmega_cState U(void)const;
   FANS3D_ThermallyPerfect_KOmega_cState U(const FANS3D_ThermallyPerfect_KOmega_pState &W);  
   
   //Inviscid Fluxes /
   FANS3D_ThermallyPerfect_KOmega_cState F(void) ;
   FANS3D_ThermallyPerfect_KOmega_cState F(void) const ;
   FANS3D_ThermallyPerfect_KOmega_cState F(const FANS3D_ThermallyPerfect_KOmega_pState &W);

  //Viscous Fluxes x y z directions/
   FANS3D_ThermallyPerfect_KOmega_cState Fv(
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   FANS3D_ThermallyPerfect_KOmega_cState Gv(
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   FANS3D_ThermallyPerfect_KOmega_cState Hv(
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
  
   static FANS3D_ThermallyPerfect_KOmega_cState FluxViscous_n(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
      const FANS3D_ThermallyPerfect_KOmega_pState &W1 ,
      const FANS3D_ThermallyPerfect_KOmega_pState &W2,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx1,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy1,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz1,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx2,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy2,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz2,
      const Vector3D &norm, const Vector3D &ts, 
      const double &deltad, const double &Volume, 
      const double &Volume_Neigbor);
   
   static FANS3D_ThermallyPerfect_KOmega_cState FluxRoe_x(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl, 
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr);
   static FANS3D_ThermallyPerfect_KOmega_cState FluxRoe_n(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
      const Vector3D &norm_dir);
   
   friend FANS3D_ThermallyPerfect_KOmega_pState HartenFixNeg(
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_a,
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_l,
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_r);
   friend FANS3D_ThermallyPerfect_KOmega_pState HartenFixPos(
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_a,
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_l,
      const FANS3D_ThermallyPerfect_KOmega_pState  &lambda_r);
   /* Eigenvalue(s), Eigenvectors (x-direction). */
   FANS3D_ThermallyPerfect_KOmega_pState lambda_x(void) ;
   FANS3D_ThermallyPerfect_KOmega_cState rc_x(const int &index) ;
   FANS3D_ThermallyPerfect_KOmega_pState lp_x(const int &index) ;
   FANS3D_ThermallyPerfect_KOmega_pState lambda_x(void) const;
   FANS3D_ThermallyPerfect_KOmega_cState rc_x(const int &index) const;
   FANS3D_ThermallyPerfect_KOmega_pState lp_x(const int &index) const;

  /* Binary arithmetic operators. */
   FANS3D_ThermallyPerfect_KOmega_pState operator +(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;
   FANS3D_ThermallyPerfect_KOmega_pState operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;
   FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a) const;
   friend FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a, const FANS3D_ThermallyPerfect_KOmega_pState &W);
   FANS3D_ThermallyPerfect_KOmega_pState operator /(const double &a) const;
   double operator *(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;
   FANS3D_ThermallyPerfect_KOmega_pState operator ^(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;
   
   /* Assignment Operator. */
  FANS3D_ThermallyPerfect_KOmega_pState& operator =(const FANS3D_ThermallyPerfect_KOmega_pState &W);

  /* Shortcut arithmetic operators. */
  FANS3D_ThermallyPerfect_KOmega_pState& operator +=(const FANS3D_ThermallyPerfect_KOmega_pState &W);
  FANS3D_ThermallyPerfect_KOmega_pState& operator -=(const FANS3D_ThermallyPerfect_KOmega_pState &W);
  
  /* Relational operators. */
  friend int operator ==(const FANS3D_ThermallyPerfect_KOmega_pState &W1,
                         const FANS3D_ThermallyPerfect_KOmega_pState &W2);
  friend int operator !=(const FANS3D_ThermallyPerfect_KOmega_pState &W1,
                         const FANS3D_ThermallyPerfect_KOmega_pState &W2);

   static FANS3D_ThermallyPerfect_KOmega_pState RoeAverage(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr);
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_x(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr);
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_x(
      const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
      const FANS3D_ThermallyPerfect_KOmega_cState &Ur);
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_n(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
      const Vector3D &norm_dir);
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_n(
      const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
      const FANS3D_ThermallyPerfect_KOmega_cState &Ur,
      const Vector3D &norm_dir);
   static FANS3D_ThermallyPerfect_KOmega_pState Reflect(
      const FANS3D_ThermallyPerfect_KOmega_pState &W,
      const Vector3D &norm_dir);
   static FANS3D_ThermallyPerfect_KOmega_pState Moving_Wall(
      const FANS3D_ThermallyPerfect_KOmega_pState &Win,
      const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
      const Vector3D &norm_dir,				 
      const Vector3D &wall_velocity,
      const Vector3D &pressure_gradient,
      const int &TEMPERATURE_BC_FLAG);
   static FANS3D_ThermallyPerfect_KOmega_pState No_Slip(
      const FANS3D_ThermallyPerfect_KOmega_pState &Win, 
      const FANS3D_ThermallyPerfect_KOmega_pState &Wout, 
      const Vector3D &norm_dir,  const Vector3D &pressure_gradient,
      const int &TEMPERATURE_BC_FLAG);
   
   // Reynolds stress tensor 
   Tensor3D lambda(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
   Tensor3D lambda(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   // Heat flux vector 
   Vector3D qflux_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
   Vector3D qflux_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const;
   
   Vector3D thermal_diffusion_t(void) const;

   static FANS3D_ThermallyPerfect_KOmega_cState Src_t(
      const FANS3D_ThermallyPerfect_KOmega_pState &Wc,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
   
   /* Index operators */
   double &operator[](int index);
   const double &operator[](int index) const;
   
   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file,
                                const FANS3D_ThermallyPerfect_KOmega_pState &W);
   friend istream& operator >> (istream &in_file,
                                FANS3D_ThermallyPerfect_KOmega_pState &W);
  
   /*  // Destructors */
/*   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; */
/*   specdata = NULL;} */
  
/*   void Deallocate(void){ if(spec != NULL) delete[] spec; */
/*   spec = NULL;  } */
  
/*   ~FANS3D_ThermallyPerfect_KOmega_pState(){ Deallocate(); } */
};


/***********************************************************
 *             FANS3D_TheramllyPerfect_cState      *
************************************************************/
class FANS3D_ThermallyPerfect_KOmega_cState : public NavierStokes3D_ThermallyPerfect_cState {
  public:

   double rhok;
   double rhoomega;
   
   Turbulence_Model_k_omega k_omega_model;

   
// constructors
   
   FANS3D_ThermallyPerfect_KOmega_cState(): NavierStokes3D_ThermallyPerfect_cState(){
      rhok = rho*ONE; rhoomega = rho*MILLION;}
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &value): NavierStokes3D_ThermallyPerfect_cState(value){
      rhok = value; rhoomega = value;}
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const double &vx,
                                  const double &vy, const double &vz,
                                  const double &En, const double &dk,
                                  const double &domega):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En ){
      rhok = dk; rhoomega = domega;}
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const double &vx,
                                  const double &vy, const double &vz,
                                  const double &En, const double &dk,
                                  const double &domega,
                                  Species *rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){
      rhok = dk; rhoomega = domega;  }
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const Vector3D &V,
                                  const double &En, const double &dk,
                                  const double &domega ):
      NavierStokes3D_ThermallyPerfect_cState(d, V, En){
      rhok = dk; rhoomega = domega;}
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const Vector3D &V,
                                  const double &En, const double &dk,
                                  const double &domega, const double &frac):
      NavierStokes3D_ThermallyPerfect_cState(d, V, En, frac){
      rhok = dk; rhoomega = domega;}
   

   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const double &vx, 
                                  const double &vy, const double &vz, 
                                  const double &En, const double &dk,
                                  const double &domega,
                                  const double &rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){
      rhok = dk; rhoomega = domega; }
   
   
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, const Vector3D &dV,
                                         const double &En,
                                         const double &dk,
                                         const double &domega,
                                         Species *rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, dV, En, rhomfrac){
      rhok = dk; rhoomega = domega;}
 //this is needed for the operator overload returns!!!!
    FANS3D_ThermallyPerfect_KOmega_cState(const FANS3D_ThermallyPerfect_KOmega_cState &U)
                { rhospec = NULL; rho = DENSITY_STDATM; set_initial_values();
                           Copy(U); }
                           

   
   FANS3D_ThermallyPerfect_KOmega_pState W(void) ; 
   FANS3D_ThermallyPerfect_KOmega_pState W(void) const;
   FANS3D_ThermallyPerfect_KOmega_pState W(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;


   /* Index operators */
   double &operator[](int index);
   const double &operator[](int index) const;
  
   //Read in ns species data
   static void set_species_data(FANS3D_ThermallyPerfect_KOmega_cState &Uo, 
                                const int &, 
                                const string *,
                                const char *, 
                                const int &, 
                                const double&,
                                const double *);

   /* VACUUM **/
   void Vacuum(){
      Euler3D_ThermallyPerfect_cState::Vacuum();
      
      rhok = ZERO;
      rhoomega = ZERO;
      
   }  
   void Copy(const FANS3D_ThermallyPerfect_KOmega_cState &U);
   double T(void) const;
   double a(void) const;     
   Vector3D thermal_diffusion_t(const double &Temp) const;
   
   double eddy_viscosity(void) const;      
   double Pr_t(void) const;      
   double Sc_t(void) const;      
   double kappa_t(void) const;     
   double Dm_t(void) const;  //molecular diffusivity
   double p(void) const;      //pressure
   double k(void) const;
   double omega(void) const;
   /* Binary arithmetic operators. */
   FANS3D_ThermallyPerfect_KOmega_cState operator +(
      const FANS3D_ThermallyPerfect_KOmega_cState &U) const;
   FANS3D_ThermallyPerfect_KOmega_cState operator -(
      const FANS3D_ThermallyPerfect_KOmega_cState &U) const;
   FANS3D_ThermallyPerfect_KOmega_cState operator *(const double &a) const;
   friend FANS3D_ThermallyPerfect_KOmega_cState operator *(
      const double &a, const FANS3D_ThermallyPerfect_KOmega_cState &U);
   FANS3D_ThermallyPerfect_KOmega_cState operator /(const double &a) const;
   
   double operator *(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;
   FANS3D_ThermallyPerfect_KOmega_cState operator ^(
      const FANS3D_ThermallyPerfect_KOmega_cState &U) const;
   
  /* Assignment Operator. */
   FANS3D_ThermallyPerfect_KOmega_cState& operator =(
      const FANS3D_ThermallyPerfect_KOmega_cState &U);
  /* Shortcut arithmetic operators. */
  FANS3D_ThermallyPerfect_KOmega_cState& operator +=(
     const FANS3D_ThermallyPerfect_KOmega_cState &U);
  FANS3D_ThermallyPerfect_KOmega_cState& operator -=(
     const FANS3D_ThermallyPerfect_KOmega_cState &U);
      
  /* Unary arithmetic operators. */
  friend FANS3D_ThermallyPerfect_KOmega_cState operator -(
     const FANS3D_ThermallyPerfect_KOmega_cState &U);
  
  /* Relational operators. */
  friend int operator ==(const FANS3D_ThermallyPerfect_KOmega_cState &U1,
                         const FANS3D_ThermallyPerfect_KOmega_cState &U2);
  friend int operator !=(const FANS3D_ThermallyPerfect_KOmega_cState &U1,
                         const FANS3D_ThermallyPerfect_KOmega_cState &U2);

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file,
                               const FANS3D_ThermallyPerfect_KOmega_cState &U);
  friend istream& operator >> (istream &in_file,
                               FANS3D_ThermallyPerfect_KOmega_cState &U);
  
/*    /\* Destructors *\/ */
/*   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata;  */
/*   specdata = NULL;  */
/*   if(Schmidt != NULL) delete[] Schmidt;  */
/*   Schmidt = NULL;  */
/*   } */
/*   void Deallocate(void){ if(rhospec != NULL) delete[] rhospec;  */
/*   rhospec = NULL;  } */
  
/*   ~FANS3D_ThermallyPerfect_KOmega_cState(){ Deallocate(); } */
  int Unphysical_Properties_Check(
     FANS3D_ThermallyPerfect_KOmega_cState &Uw,
     FANS3D_ThermallyPerfect_KOmega_cState &Ue,
     FANS3D_ThermallyPerfect_KOmega_cState &Ut,
     FANS3D_ThermallyPerfect_KOmega_cState &Ub,
     FANS3D_ThermallyPerfect_KOmega_cState &Un,
     FANS3D_ThermallyPerfect_KOmega_cState &Us, 
     const int n);
  
};

inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::U(void){
  
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   Temp.rhoomega = rho*omega;


   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
  
}

inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::U(void)const{
  
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   Temp.rhoomega = rho*omega;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
  
}

inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::U(const FANS3D_ThermallyPerfect_KOmega_pState &W){
  if(ns == W.ns){ //check that species are equal
    FANS3D_ThermallyPerfect_KOmega_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    Temp.rhok = W.rho*W.k;
    Temp.rhoomega = W.rho*W.omega;
    for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
    }
 
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

inline FANS3D_ThermallyPerfect_KOmega_cState U(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   Temp.rhok = W.rho*W.k;
   Temp.rhoomega = W.rho*W.omega;
   for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
   }
   
   return Temp;
}



inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::W(void){
  
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();
   Temp.omega = omega();
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
   
   
}
//Copy
inline void FANS3D_ThermallyPerfect_KOmega_pState::Copy(
   const FANS3D_ThermallyPerfect_KOmega_pState &W){
   rho = W.rho;
   v = W.v; 
   p = W.p;  
   k = W.k;
   omega = W.omega;

   for( int i=0; i<ns; i++){
      spec[i] = W.spec[i];
   }
}

//----------------- Index Operator ------------------------/
inline double& FANS3D_ThermallyPerfect_KOmega_pState::operator[](int index) {
   
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   case 6:
      return k;
   case 7:
      return omega;

   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;
   };
}

inline const double& FANS3D_ThermallyPerfect_KOmega_pState::operator[](int index) const {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   case 6:
      return k;
   case 7:
      return omega;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;     
      
   };


}

inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::W(void)const{
  
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();  
   Temp.omega = omega();
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
   
   
}

inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::W(
   const FANS3D_ThermallyPerfect_KOmega_cState &U) const{
  if(ns == U.ns){ //check that species are equal   
    FANS3D_ThermallyPerfect_KOmega_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    Temp.p = U.p();
    Temp.k = U.k();  
    Temp.omega = U.omega();
  
    for(int i=0; i<U.ns; i++){
      Temp.spec[i] = U.rhospec[i]/U.rho;
      Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
      Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
    }
  
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  } 
}

inline FANS3D_ThermallyPerfect_KOmega_pState W(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
  FANS3D_ThermallyPerfect_KOmega_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  Temp.p = U.p();
  Temp.k = U.k();
  Temp.omega = U.omega();
     
  for(int i=0; i<U.ns; i++){
    Temp.spec[i] = U.rhospec[i]/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  }
 
  return Temp;
}

// /******************************************************
//  Calculating the thermal diffusion component of 
//  the heat flux vector (qflux)
//   sum( hs * Ds * grad cs)
// *******************************************************/

inline Vector3D FANS3D_ThermallyPerfect_KOmega_pState::thermal_diffusion_t(void) const{
   Vector3D sum;
   sum.zero();
   double Temp = T();
   //problems with Species overloaded operators
   for(int i=0; i<ns; i++){
      sum  +=  (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
         * rho*Dm_t() * spec[i].gradc;
   }

   
   return sum;
}

inline Vector3D FANS3D_ThermallyPerfect_KOmega_cState::thermal_diffusion_t(
   const double &Temp) const{
   Vector3D sum;
   sum.zero();
   
   for(int i=0; i<ns; i++){
      sum  +=   (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
         * rho*Dm_t()*rhospec[i].gradc;
      
   }
   return sum/(rho*rho);
   
}

/* Copy */
inline void FANS3D_ThermallyPerfect_KOmega_cState::Copy(
   const FANS3D_ThermallyPerfect_KOmega_cState &U){
   rho = U.rho;
   rhov = U.rhov; 
   E = U.E; 
   rhok = U.rhok;
   rhoomega = U.rhoomega;
   
   for( int i=0; i<ns; i++){ 
      rhospec[i] = U.rhospec[i];
   } 
   
}
//index operators
inline double& FANS3D_ThermallyPerfect_KOmega_cState::operator[](int index) {

   switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   case 6:
      return rhok;
   case 7:
      return rhoomega;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;     
      
  
   };
    
}

inline const double& FANS3D_ThermallyPerfect_KOmega_cState::operator[](int index) const{
  
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   case 6:
      return rhok;
   case 7:
      return rhoomega;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;     
 
   };
 
}
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED 

/****************** NAVIERSTOKES3D_ThermallyPerfectState.h *************
  This class defines the state variables and constructors for the 
  NavierStokes3D thermally perfect gaseous mixture class.
 ***********************************************************************/

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED 
#define _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

class NavierStokes3D_ThermallyPerfect_cState;
class NavierStokes3D_ThermallyPerfect_pState;

// Required C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

// Include required CFFC header files

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../Euler/Euler3DThermallyPerfectState.h"
#endif  //EULER3D_THERMALLYPERFECT_STATE_INCLUDE

#ifndef _REACTIONS_INCLUDED
#include "../Reactions/Reactions.h"
#endif // _REACTIONS_INCLUDED

/********************************************************
 * Class: NavierStokes3D_ThermallyPerfect_pState        *
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
  NASARP1311data *specdata:   Global Species Data

  num_vars:         number of total variables (5+ns)
  Reaction_set React:         Global Reaction Data

  low_temp_range:      low temp data range
  high_temp_range:     high temp data range
***********************************************************/

class NavierStokes3D_ThermallyPerfect_pState : public Euler3D_ThermallyPerfect_pState {
  public:
   // constructor
   NavierStokes3D_ThermallyPerfect_pState(): Euler3D_ThermallyPerfect_pState(){ }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &value):
      Euler3D_ThermallyPerfect_pState(value){ }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const Vector3D &V,
                                          const double &pre) :
      Euler3D_ThermallyPerfect_pState(d, V, pre){  }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const Vector3D &V,
                                          const double &pre, const double &frac):
      Euler3D_ThermallyPerfect_pState(d, V, pre, frac){ }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const double &vx,
                                          const double &vy, const double &vz,
                                          const double &pre, const double &frac):
      Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac){ }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const double &vx,
                                          const double &vy, const double &vz,
                                          const double &pre):
      Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre) { }
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const double &vx,
                                          const double &vy,  const double &vz,
                                          const double &pre, Species *mfrac):
      Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac){ }
     
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const Vector3D &V,
                                          const double &pre, Species *mfrac):
      Euler3D_ThermallyPerfect_pState(d, V, pre, mfrac){  }
   
   
   NavierStokes3D_ThermallyPerfect_pState(const double &d, const Vector3D &V,
                                          const double &pre, double *mfrac):
      Euler3D_ThermallyPerfect_pState(d, V, pre, mfrac){
      set_initial_values(mfrac); }

    //this is needed for the operator overload returns!!!!
   NavierStokes3D_ThermallyPerfect_pState(const  NavierStokes3D_ThermallyPerfect_pState &W)
                {spec = NULL; rho = DENSITY_STDATM; set_initial_values(); Copy(W);}
                

   // Conserved solution state. /
   NavierStokes3D_ThermallyPerfect_cState U(void);
   NavierStokes3D_ThermallyPerfect_cState U(void)const;
   NavierStokes3D_ThermallyPerfect_cState U(const NavierStokes3D_ThermallyPerfect_pState &W);  
    
  //Inviscid Fluxes /
   NavierStokes3D_ThermallyPerfect_cState F(void) ;
   NavierStokes3D_ThermallyPerfect_cState F(void) const ;
   NavierStokes3D_ThermallyPerfect_cState F(const NavierStokes3D_ThermallyPerfect_pState &W);
 
  //Viscous Fluxes x y z directions/
  NavierStokes3D_ThermallyPerfect_cState Fv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) ;
  NavierStokes3D_ThermallyPerfect_cState Gv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) ;
  
  NavierStokes3D_ThermallyPerfect_cState Hv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) ;
  NavierStokes3D_ThermallyPerfect_cState Fv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
  NavierStokes3D_ThermallyPerfect_cState Gv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
  NavierStokes3D_ThermallyPerfect_cState Hv(
     const NavierStokes3D_ThermallyPerfect_pState &dWdx,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
  
  static NavierStokes3D_ThermallyPerfect_cState FluxViscous_n(
     const NavierStokes3D_ThermallyPerfect_pState &Wl,
     const NavierStokes3D_ThermallyPerfect_pState &Wr,
     const NavierStokes3D_ThermallyPerfect_pState &W1 ,
     const NavierStokes3D_ThermallyPerfect_pState &W2,
     const NavierStokes3D_ThermallyPerfect_pState &dWdx1,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy1,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz1,
     const NavierStokes3D_ThermallyPerfect_pState &dWdx2,
     const NavierStokes3D_ThermallyPerfect_pState &dWdy2,
     const NavierStokes3D_ThermallyPerfect_pState &dWdz2,
     const Vector3D &norm, const Vector3D &ts, 
     const double &deltad, const double &Volume, const double &Volume_Neigbor);
  
  static NavierStokes3D_ThermallyPerfect_cState FluxRoe_x(
     const NavierStokes3D_ThermallyPerfect_pState &Wl, 
     const NavierStokes3D_ThermallyPerfect_pState &Wr);
  static NavierStokes3D_ThermallyPerfect_cState FluxRoe_n(
     const NavierStokes3D_ThermallyPerfect_pState &Wl,
     const NavierStokes3D_ThermallyPerfect_pState &Wr,
     const Vector3D &norm_dir);

  friend NavierStokes3D_ThermallyPerfect_pState HartenFixNeg(
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_a,
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_l,
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_r);
  friend NavierStokes3D_ThermallyPerfect_pState HartenFixPos(
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_a,
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_l,
     const NavierStokes3D_ThermallyPerfect_pState  &lambda_r);
  /* Eigenvalue(s), Eigenvectors (x-direction). */
  NavierStokes3D_ThermallyPerfect_pState lambda_x(void) ;
  NavierStokes3D_ThermallyPerfect_cState rc_x(const int &index) ;
  NavierStokes3D_ThermallyPerfect_pState lp_x(const int &index) ;
  NavierStokes3D_ThermallyPerfect_pState lambda_x(void) const;
  NavierStokes3D_ThermallyPerfect_cState rc_x(const int &index) const;
  NavierStokes3D_ThermallyPerfect_pState lp_x(const int &index) const;

  /* Binary arithmetic operators. */
  NavierStokes3D_ThermallyPerfect_pState operator +(const NavierStokes3D_ThermallyPerfect_pState &W) const;
  NavierStokes3D_ThermallyPerfect_pState operator -(const NavierStokes3D_ThermallyPerfect_pState &W) const;
  NavierStokes3D_ThermallyPerfect_pState operator *(const double &a) const;
  friend NavierStokes3D_ThermallyPerfect_pState operator *(const double &a, const NavierStokes3D_ThermallyPerfect_pState &W);
  NavierStokes3D_ThermallyPerfect_pState operator /(const double &a) const;
  double operator *(const NavierStokes3D_ThermallyPerfect_pState &W) const;
  NavierStokes3D_ThermallyPerfect_pState operator ^(const NavierStokes3D_ThermallyPerfect_pState &W) const;

  /* Assignment Operator. */
  NavierStokes3D_ThermallyPerfect_pState& operator =(const NavierStokes3D_ThermallyPerfect_pState &W);

  /* Shortcut arithmetic operators. */
  NavierStokes3D_ThermallyPerfect_pState& operator +=(const NavierStokes3D_ThermallyPerfect_pState &W);
  NavierStokes3D_ThermallyPerfect_pState& operator -=(const NavierStokes3D_ThermallyPerfect_pState &W);
  
  /* Relational operators. */
  friend int operator ==(const NavierStokes3D_ThermallyPerfect_pState &W1,
                         const NavierStokes3D_ThermallyPerfect_pState &W2);
  friend int operator !=(const NavierStokes3D_ThermallyPerfect_pState &W1,
                         const NavierStokes3D_ThermallyPerfect_pState &W2);

  static NavierStokes3D_ThermallyPerfect_pState RoeAverage(
     const NavierStokes3D_ThermallyPerfect_pState &Wl,
     const NavierStokes3D_ThermallyPerfect_pState &Wr);
  static NavierStokes3D_ThermallyPerfect_cState FluxHLLE_x(
      const NavierStokes3D_ThermallyPerfect_pState &Wl,
      const NavierStokes3D_ThermallyPerfect_pState &Wr);
  static NavierStokes3D_ThermallyPerfect_cState FluxHLLE_x(
     const NavierStokes3D_ThermallyPerfect_cState &Ul,
     const NavierStokes3D_ThermallyPerfect_cState &Ur);
  static NavierStokes3D_ThermallyPerfect_cState FluxHLLE_n(
      const NavierStokes3D_ThermallyPerfect_pState &Wl,
      const NavierStokes3D_ThermallyPerfect_pState &Wr,
      const Vector3D &norm_dir);
  static NavierStokes3D_ThermallyPerfect_cState FluxHLLE_n(
      const NavierStokes3D_ThermallyPerfect_cState &Ul,
      const NavierStokes3D_ThermallyPerfect_cState &Ur,
      const Vector3D &norm_dir);
  static NavierStokes3D_ThermallyPerfect_pState Reflect(
      const NavierStokes3D_ThermallyPerfect_pState &W,
      const Vector3D &norm_dir);
  static NavierStokes3D_ThermallyPerfect_pState Moving_Wall(
     const NavierStokes3D_ThermallyPerfect_pState &Win,
     const NavierStokes3D_ThermallyPerfect_pState &Wout,
     const Vector3D &norm_dir,				 
     const Vector3D &wall_velocity,
     const Vector3D &pressure_gradient,
     const int &TEMPERATURE_BC_FLAG);
  static NavierStokes3D_ThermallyPerfect_pState No_Slip(
     const NavierStokes3D_ThermallyPerfect_pState &Win, 
     const NavierStokes3D_ThermallyPerfect_pState &Wout,
     const Vector3D &norm_dir, 
     const Vector3D &pressure_gradient,
     const int &TEMPERATURE_BC_FLAG);
   
   // shear stress tensor 
   Tensor3D tau(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   Tensor3D tau(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
   // Heat flux vector 
   Vector3D qflux(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   Vector3D qflux(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
   
   Vector3D thermal_diffusion(void) const;
   /* Input-output operators. */
 

   /****** Source terms associated with finite-rate chemistry ******/
   NavierStokes3D_ThermallyPerfect_cState Sw(int &REACT_SET_FLAG, const int flow_type) const;
   
   friend ostream& operator << (ostream &out_file,
                                const NavierStokes3D_ThermallyPerfect_pState &W);
   friend istream& operator >> (istream &in_file,
                                NavierStokes3D_ThermallyPerfect_pState &W);
  
   //Read in ns species data
/*    void set_species_data(const int &n, */
/*                          const string *S, */
/*                          const char *PATH, */
/*                          const int &debug,  */
/*                          const double &Mr,  */
/*                          const double* Sc, */
/*                          const int &trans_data); */

};


/***********************************************************
 *             NavierStokes3D_TheramllyPerfect_cState      *
************************************************************/
class NavierStokes3D_ThermallyPerfect_cState : public Euler3D_ThermallyPerfect_cState {
  public:

// constructors
   
   NavierStokes3D_ThermallyPerfect_cState(): Euler3D_ThermallyPerfect_cState(){ }
 
   NavierStokes3D_ThermallyPerfect_cState(const double &value): Euler3D_ThermallyPerfect_cState(value){ }
  
   NavierStokes3D_ThermallyPerfect_cState(const double &d, const double &vx,
                                          const double &vy, const double &vz,
                                          const double &En ):
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En ){ }
   
   NavierStokes3D_ThermallyPerfect_cState(const double &d, const double &vx,
                                          const double &vy, const double &vz,
                                          const double &En,
                                          Species *rhomfrac):
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){ }
   
   NavierStokes3D_ThermallyPerfect_cState(const double &d, const Vector3D &V,
                                          const double &En):
      Euler3D_ThermallyPerfect_cState(d, V, En){ }
   
   NavierStokes3D_ThermallyPerfect_cState(const double &d, const Vector3D &V,
                                          const double &En, const double &frac):
      Euler3D_ThermallyPerfect_cState(d, V, En, frac){}
   

   NavierStokes3D_ThermallyPerfect_cState(const double &d, const double &vx, 
                                          const double &vy, const double &vz, 
                                          const double &En,	
                                          const double &rhomfrac):
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){  }
   
   
   NavierStokes3D_ThermallyPerfect_cState(const double &d, const Vector3D &dV,
                                          const double &En, Species *rhomfrac):
      Euler3D_ThermallyPerfect_cState(d, dV, En, rhomfrac){ }

 //this is needed for the operator overload returns!!!!
       NavierStokes3D_ThermallyPerfect_cState(const NavierStokes3D_ThermallyPerfect_cState &U)
               { rhospec = NULL; rho = DENSITY_STDATM; set_initial_values();
                          Copy(U); }
                          

   
   NavierStokes3D_ThermallyPerfect_pState W(void) ; 
   NavierStokes3D_ThermallyPerfect_pState W(void) const;
   NavierStokes3D_ThermallyPerfect_pState W(const NavierStokes3D_ThermallyPerfect_cState &U) const;
   friend NavierStokes3D_ThermallyPerfect_pState W(const NavierStokes3D_ThermallyPerfect_cState &U); 

/*   /\* Index operators *\/ */
/*   double &operator[](int index); */
/*   const double &operator[](int index) const; */
 
  /* Binary arithmetic operators. */
  NavierStokes3D_ThermallyPerfect_cState operator +(
     const NavierStokes3D_ThermallyPerfect_cState &U) const;
  NavierStokes3D_ThermallyPerfect_cState operator -(
     const NavierStokes3D_ThermallyPerfect_cState &U) const;
  NavierStokes3D_ThermallyPerfect_cState operator *(const double &a) const;
  friend NavierStokes3D_ThermallyPerfect_cState operator *(
     const double &a, const NavierStokes3D_ThermallyPerfect_cState &U);
  NavierStokes3D_ThermallyPerfect_cState operator /(const double &a) const;
  
  double operator *(const NavierStokes3D_ThermallyPerfect_cState &U) const;
  NavierStokes3D_ThermallyPerfect_cState operator ^(
     const NavierStokes3D_ThermallyPerfect_cState &U) const;

  /* Assignment Operator. */
   NavierStokes3D_ThermallyPerfect_cState& operator =(
      const NavierStokes3D_ThermallyPerfect_cState &U);
  /* Shortcut arithmetic operators. */
  NavierStokes3D_ThermallyPerfect_cState& operator +=(
     const NavierStokes3D_ThermallyPerfect_cState &U);
  NavierStokes3D_ThermallyPerfect_cState& operator -=(
     const NavierStokes3D_ThermallyPerfect_cState &U);
      
  /* Unary arithmetic operators. */
  friend NavierStokes3D_ThermallyPerfect_cState operator -(
     const NavierStokes3D_ThermallyPerfect_cState &U);
  
  /* Relational operators. */
  friend int operator ==(const NavierStokes3D_ThermallyPerfect_cState &U1,
                         const NavierStokes3D_ThermallyPerfect_cState &U2);
  friend int operator !=(const NavierStokes3D_ThermallyPerfect_cState &U1,
                         const NavierStokes3D_ThermallyPerfect_cState &U2);

  Vector3D thermal_diffusion(const double &T) const;
  
  //Read in ns species data
/*   void set_species_data(const int &n, */
/*                         const string *S, */
/*                         const char *PATH, */
/*                         const int &debug,  */
/*                         const double &Mr,  */
/*                         const double* Sc, */
/*                         const int &trans_data); */
  
};

// Can not use ones from NavierStokes3D _ThermallyPerfect_pState Class due to the return type  
inline NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::U(void){
  
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
  
}


inline NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::U(void)const{
  
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
  
}

inline NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::U(const NavierStokes3D_ThermallyPerfect_pState &W){
  if(ns == W.ns){ //check that species are equal
    NavierStokes3D_ThermallyPerfect_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
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

inline NavierStokes3D_ThermallyPerfect_cState U(const NavierStokes3D_ThermallyPerfect_pState &W) {
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
   }
   
   return Temp;
}



inline NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_cState::W(void){
  
   NavierStokes3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
   
   
}


inline NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_cState::W(void)const{
  
   NavierStokes3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
   
   
}

inline NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_cState::W(
   const NavierStokes3D_ThermallyPerfect_cState &U) const{
  if(ns == U.ns){ //check that species are equal   
    NavierStokes3D_ThermallyPerfect_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    Temp.p = U.p();
  
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

inline NavierStokes3D_ThermallyPerfect_pState W(const NavierStokes3D_ThermallyPerfect_cState &U) {
  NavierStokes3D_ThermallyPerfect_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  Temp.p = U.p();

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

inline Vector3D NavierStokes3D_ThermallyPerfect_pState::thermal_diffusion(void) const{
   Vector3D sum;
   sum.zero();
   double Temp = T();
   //problems with Species overloaded operators
   for(int i=0; i<ns; i++){
      sum  +=  (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
         * spec[i].diffusion_coef * spec[i].gradc;
   }
   return sum;
}


inline Vector3D NavierStokes3D_ThermallyPerfect_cState::thermal_diffusion(
   const double &Temp) const{
   Vector3D sum;
   sum.zero();
   
   for(int i=0; i<ns; i++){
      sum  +=   (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
         * rhospec[i].diffusion_coef*rhospec[i].gradc;
      
   }
   return sum/(rho*rho);
   
}

#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED 

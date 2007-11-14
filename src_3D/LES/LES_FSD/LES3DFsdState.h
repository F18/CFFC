/******************** LES3DFsdState.h *************************
  This class defines the state variables and constructors 
  for LES Euler3d thermally perfect mixtures
***************************************************************/
#ifndef _LES3DFSD_STATE_INCLUDED 
#define _LES3DFSD_STATE_INCLUDED

class LES3DFsd_cState;
class LES3DFsd_pState;

// Required C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

// Required CFFC header files
#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../../NavierStokes/NavierStokes3DThermallyPerfectState.h"
#endif  //_NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

#ifndef  _TURBULENCE_MODELLING_INCLUDED 
#include "../../TurbulenceModelling/TurbulenceModelling.h"
#endif  //_TURBULENCE_MODELLING_INCLUDED

#define NUM_LES3D_VAR_EXTRA 3  //k,C,Fsd

/********************************************************
 * Class: LES3DFsd_pState                               *
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

  NUM_VAR_3D;         number of total variables (5+nextra+ns)
  Reaction_set React;         Global Reaction Data

  low_temp_range;      low temp data range
  high_temp_range;     high temp data range

***********************************************************/
class LES3DFsd_pState :public NavierStokes3D_ThermallyPerfect_pState {
  public:
   double C;
   double Fsd;
   double k;

   // constructor
   LES3DFsd_pState(): 
     NavierStokes3D_ThermallyPerfect_pState(){C = ONE; Fsd = MILLION; k = ONE;}
   
   LES3DFsd_pState(const double &value):
     NavierStokes3D_ThermallyPerfect_pState(value){C = value; Fsd = value; k = value;}
   
   LES3DFsd_pState(const double &d, const Vector3D &V,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre){C = CC; Fsd = FFsd; k = K;}

   LES3DFsd_pState(const double &d, const Vector3D &V,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K,
                   const double &frac):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, frac){C = CC; Fsd = FFsd; k = K;}
   
   LES3DFsd_pState(const double &d, const double &vx,
                   const double &vy, const double &vz,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K,
                   const double &frac):
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac){C = CC; Fsd = FFsd; k = K;}
   
   LES3DFsd_pState(const double &d, const double &vx,
                   const double &vy, const double &vz,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K):
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre) {C = CC; Fsd = FFsd; k = K;}
   
   LES3DFsd_pState(const double &d, const double &vx,
                   const double &vy, const double &vz,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K,
                   Species *mfrac):
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac){C = CC; Fsd = FFsd; k = K;}
   
   LES3DFsd_pState(const double &d, const Vector3D &V,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K,
                   Species *mfrac):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){C = CC; Fsd = FFsd; k = K;}
   
   LES3DFsd_pState(const double &d, const Vector3D &V,
                   const double &pre, double *mfrac):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){C = ONE; Fsd = MILLION; k = ONE; set_initial_values(mfrac);}
      
   LES3DFsd_pState(const double &d, const Vector3D &V,
                   const double &pre, const double &CC, 
                   const double &FFsd, const double &K,
                   double *mfrac):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac){C = CC; Fsd = FFsd; k = K; set_initial_values(mfrac);}
  
// this is needed for the operator overload returns!!!!
  LES3DFsd_pState(const  LES3DFsd_pState &W) {spec = NULL; rho = DENSITY_STDATM; set_initial_values(); Copy(W);}

   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);

   void Vacuum(){ Euler3D_ThermallyPerfect_pState::Vacuum(); C = ZERO; Fsd = ZERO; k = ZERO;}

   void Copy(const LES3DFsd_pState &W);

   double E(void) const ;      //mixture total internal energy
   double H(void) const ;      //mixture total enthalpy
   double Hs(void) const;
   double pmodified(void) const;
   double a(void);                
   double a(void) const;
   static double Mref;               // Mref for Precondtioning (normally set to incoming freestream Mach) 
/*    double Laminar_Flame_Speed; */
/*    double Laminar_Flame_Thickness; */
/*    double Adiabatic_Temperature; */
/*    double Equivalence_Ratio; */
/*    double Reactants_Density; */
/*****************************************************************/
   Tensor3D rotation_tensor(const LES3DFsd_pState &dWdx, 
                            const LES3DFsd_pState &dWdy, 
                            const LES3DFsd_pState &dWdz) const;
   double abs_strain_rate(const LES3DFsd_pState &dWdx, 
                          const LES3DFsd_pState &dWdy, 
                          const LES3DFsd_pState &dWdz) const;
   double eddy_viscosity(const LES3DFsd_pState &dWdx, 
                         const LES3DFsd_pState &dWdy, 
                         const LES3DFsd_pState &dWdz,
                         const int &Flow_Type, 
                         const double &Volume) const;      
   double Pr_t(void) const;      
   double Sc_t(void) const;      
   Vector3D thermal_diffusion_t(const LES3DFsd_pState &dWdx,
 			        const LES3DFsd_pState &dWdy,
			        const LES3DFsd_pState &dWdz,
                                const int &Flow_Type, 
                                const double &Volume) const;
/************** Premixed combustion ****************************/
   LES3DFsd_pState premixed_mfrac(const LES3DFsd_pState &Wo, const LES3DFsd_pState &W);
/************** C-FSD-k Model Source Term **********************/
   double filter_width(const double &Volume) const;
   double HeatRelease_Parameter(void)const;
   double SFS_Kinetic_Energy_Fsd(const LES3DFsd_pState &dWdx,
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const int &Flow_Type,
                                 const double &Volume) const;
   double Efficiency_Function_Fsd(const LES3DFsd_pState &dWdx,
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
                                  const int &Flow_Type,
                                  const double &Volume) const; 
   double Progvar_Species_Grad(void) const;
   double Reaction_Rate_Fsd(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const;
   double M_x(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;
   double M_y(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;
   double M_z(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;
   double Resolved_Strain(const LES3DFsd_pState &dWdx,
                          const LES3DFsd_pState &dWdy,
                          const LES3DFsd_pState &dWdz) const;
   double Resolved_Propagation_Curvature(const LES3DFsd_pState &dWdx,
                                         const LES3DFsd_pState &dWdy,
                                         const LES3DFsd_pState &dWdz) const;
   double SFS_Strain(const LES3DFsd_pState &dWdx,
                     const LES3DFsd_pState &dWdy,
                     const LES3DFsd_pState &dWdz,
                     const int &Flow_Type,
                     const double &Volume) const;
   double SFS_Curvature(const LES3DFsd_pState &dWdx,
                        const LES3DFsd_pState &dWdy,
                        const LES3DFsd_pState &dWdz) const;
   double M_xx(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdx_dx,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdx_dz) const;
   double M_yy(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdy_dy,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdy_dz) const;
   double M_zz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;
   double M_xy(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdy_dy,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdy_dz) const;
   double M_xz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;
   double M_yz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;
   double Resolved_Curvature(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdx_dx,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const;
   double Resolved_Propagation(const LES3DFsd_pState &dWdx,
                               const LES3DFsd_pState &dWdy,
                               const LES3DFsd_pState &dWdz,
                               const LES3DFsd_pState &d_dWdx_dx,
                               const LES3DFsd_pState &d_dWdy_dy,
                               const LES3DFsd_pState &d_dWdz_dz,
                               const LES3DFsd_pState &d_dWdx_dy,
                               const LES3DFsd_pState &d_dWdx_dz,
                               const LES3DFsd_pState &d_dWdy_dz) const;
   double Resolved_Convection_Progvar (const LES3DFsd_pState &dWdx,
                                       const LES3DFsd_pState &dWdy,
                                       const LES3DFsd_pState &dWdz) const;
   double Resolved_Convection_Fsd (const LES3DFsd_pState &dWdx,
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz) const;
   double NGT_Progvar (const LES3DFsd_pState &dWdx,
                       const LES3DFsd_pState &dWdy,
                       const LES3DFsd_pState &dWdz) const;
   double NGT_Fsd (const LES3DFsd_pState &dWdx,
                   const LES3DFsd_pState &dWdy,
                   const LES3DFsd_pState &dWdz,
                   const LES3DFsd_pState &d_dWdx_dx,
                   const LES3DFsd_pState &d_dWdy_dy,
                   const LES3DFsd_pState &d_dWdz_dz,
                   const LES3DFsd_pState &d_dWdx_dy,
                   const LES3DFsd_pState &d_dWdx_dz,
                   const LES3DFsd_pState &d_dWdy_dz) const;
   double SFS_Diffusion_Progvar (const LES3DFsd_pState &dWdx,
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const LES3DFsd_pState &d_dWdx_dx,
                                 const LES3DFsd_pState &d_dWdy_dy,
                                 const LES3DFsd_pState &d_dWdz_dz,
                                 const LES3DFsd_pState &d_dWdx_dy,
                                 const LES3DFsd_pState &d_dWdx_dz,
                                 const LES3DFsd_pState &d_dWdy_dz,
                                 const int &Flow_Type,
                                 const double &Volume) const;
   double SFS_Diffusion_Fsd (const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdx_dx,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdx_dz,
             	             const LES3DFsd_pState &d_dWdy_dz,
                             const int &Flow_Type,
                             const double &Volume) const; 
   double Heat_Release_Strain (const LES3DFsd_pState &dWdx,
                               const LES3DFsd_pState &dWdy,
                               const LES3DFsd_pState &dWdz,
                               const LES3DFsd_pState &d_dWdx_dx,
                               const LES3DFsd_pState &d_dWdy_dy,
                               const LES3DFsd_pState &d_dWdz_dz,
                               const LES3DFsd_pState &d_dWdx_dy,
                               const LES3DFsd_pState &d_dWdx_dz,
                               const LES3DFsd_pState &d_dWdy_dz) const;
   double Net_Rate_Change_Progvar (const LES3DFsd_pState &dWdx,
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz,
                                   const LES3DFsd_pState &d_dWdx_dx,
                                   const LES3DFsd_pState &d_dWdy_dy,
                                   const LES3DFsd_pState &d_dWdz_dz,
                                   const LES3DFsd_pState &d_dWdx_dy,
                                   const LES3DFsd_pState &d_dWdx_dz,
		                   const LES3DFsd_pState &d_dWdy_dz,
                                   const int &Flow_Type,
                                   const double &Volume) const;
   double Net_Rate_Change_Fsd (const LES3DFsd_pState &dWdx,
                               const LES3DFsd_pState &dWdy,
                               const LES3DFsd_pState &dWdz,
                               const LES3DFsd_pState &d_dWdx_dx,
                               const LES3DFsd_pState &d_dWdy_dy,
                               const LES3DFsd_pState &d_dWdz_dz,
                               const LES3DFsd_pState &d_dWdx_dy,
                               const LES3DFsd_pState &d_dWdx_dz,
		               const LES3DFsd_pState &d_dWdy_dz,
                               const int &Flow_Type,
                               const double &Volume) const;
   double K_equ_sources(const LES3DFsd_pState &dWdx,
                        const LES3DFsd_pState &dWdy,
                        const LES3DFsd_pState &dWdz,
                        const int &Flow_Type,
                        const double &Volume) const;
 
   //Conserved solution state.
   LES3DFsd_cState U(void);
   LES3DFsd_cState U(void)const;
   LES3DFsd_cState U(const LES3DFsd_pState &W);  
   
   //Inviscid Fluxes /
   LES3DFsd_cState F(void) ;
   LES3DFsd_cState F(void) const ;
   LES3DFsd_cState F(const LES3DFsd_pState &W);

  //Viscous Fluxes x y z directions
   LES3DFsd_cState Fv(const LES3DFsd_pState &dWdx,
                      const LES3DFsd_pState &dWdy,
                      const LES3DFsd_pState &dWdz,
                      const int &Flow_Type,
                      const double &Volume) const;
   LES3DFsd_cState Gv(const LES3DFsd_pState &dWdx,
                      const LES3DFsd_pState &dWdy,
                      const LES3DFsd_pState &dWdz,
                      const int &Flow_Type,
                      const double &Volume) const;
   LES3DFsd_cState Hv(const LES3DFsd_pState &dWdx,
                      const LES3DFsd_pState &dWdy,
                      const LES3DFsd_pState &dWdz,
                      const int &Flow_Type,
                      const double &Volume) const;
   static LES3DFsd_cState FluxViscous_n(const LES3DFsd_pState &Wl,
                                        const LES3DFsd_pState &Wr,
                                        const LES3DFsd_pState &W1,
                                        const LES3DFsd_pState &W2,
                                        const LES3DFsd_pState &dWdx1,
                                        const LES3DFsd_pState &dWdy1,
                                        const LES3DFsd_pState &dWdz1,
                                        const LES3DFsd_pState &dWdx2,
                                        const LES3DFsd_pState &dWdy2,
                                        const LES3DFsd_pState &dWdz2,
                                        const Vector3D &norm, const Vector3D &ts, 
                                        const double &deltad, const double &Volume, 
                                        const double &Volume_Neigbor, const int &Flow_Type);
   static LES3DFsd_cState FluxRoe_x(const LES3DFsd_pState &Wl, 
                                    const LES3DFsd_pState &Wr);
   static LES3DFsd_cState FluxRoe_n(const LES3DFsd_pState &Wl,
                                    const LES3DFsd_pState &Wr,
                                    const Vector3D &norm_dir);
   friend LES3DFsd_pState HartenFixNeg(const LES3DFsd_pState  &lambda_a,
                                       const LES3DFsd_pState  &lambda_l,
                                       const LES3DFsd_pState  &lambda_r);
   friend LES3DFsd_pState HartenFixPos(const LES3DFsd_pState  &lambda_a,
                                       const LES3DFsd_pState  &lambda_l,
                                       const LES3DFsd_pState  &lambda_r);
   // Eigenvalue(s), Eigenvectors (x-direction)
   LES3DFsd_pState lambda_x(void) ;
   LES3DFsd_cState rc_x(const int &index) ;
   LES3DFsd_pState lp_x(const int &index) ;
   LES3DFsd_pState lambda_x(void) const;
   LES3DFsd_cState rc_x(const int &index) const;
   LES3DFsd_pState lp_x(const int &index) const;

/*************** Preconditioner ****************************/
  LES3DFsd_cState rc_x_precon(const int &index,const double &MR2) const; 
  LES3DFsd_pState lp_x_precon(const int &index,const double &MR2) const; 
  double u_plus_aprecon(const double &u, double &deltax, const double &lengthx, const double &dTime) const;
  void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
  double Mr2(const double &deltax, const double &lengthx, const double &dTime) const;
  void Low_Mach_Number_Preconditioner(DenseMatrix &P, const double &deltax, const double &lengthx, const double &dTime) const; 
  void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv, const double &deltax, const double &lengthx, const double &dTime) const; 

  /* compute the dWdU */
  void dWdU(DenseMatrix &dWdU);

/************************************************************/
  void SemiImplicitSourceJacobi(const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const double &d_dWdx_dW, 
                                const double &d_dWdy_dW,
                                const double &d_dWdz_dW,
                                DenseMatrix &dStdW,
                                const int &Flow_Type,
                                const double &Volume)const;

  /* Binary arithmetic operators. */
   LES3DFsd_pState operator +(const LES3DFsd_pState &W) const;
   LES3DFsd_pState operator -(const LES3DFsd_pState &W) const;
   LES3DFsd_pState operator *(const double &a) const;
   friend LES3DFsd_pState operator *(const double &a, const LES3DFsd_pState &W);
   LES3DFsd_pState operator /(const double &a) const;
   double operator *(const LES3DFsd_pState &W) const;
   LES3DFsd_pState operator ^(const LES3DFsd_pState &W) const;
   
   /* Assignment Operator. */
  LES3DFsd_pState& operator =(const LES3DFsd_pState &W);

  /* Shortcut arithmetic operators. */
  LES3DFsd_pState& operator +=(const LES3DFsd_pState &W);
  LES3DFsd_pState& operator -=(const LES3DFsd_pState &W);
  
  /* Relational operators. */
  friend int operator ==(const LES3DFsd_pState &W1,
                         const LES3DFsd_pState &W2);
  friend int operator !=(const LES3DFsd_pState &W1,
                         const LES3DFsd_pState &W2);

   static LES3DFsd_pState RoeAverage(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr);
   static LES3DFsd_cState FluxHLLE_x(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr);
   static LES3DFsd_cState FluxHLLE_x(const LES3DFsd_cState &Ul,
                                     const LES3DFsd_cState &Ur);
   static LES3DFsd_cState FluxHLLE_n(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr,
                                     const Vector3D &norm_dir);
   static LES3DFsd_cState FluxHLLE_n(const LES3DFsd_cState &Ul,
                                     const LES3DFsd_cState &Ur,
                                     const Vector3D &norm_dir);

/****************************************************************/
   static LES3DFsd_cState FluxAUSMplus_up_x(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr);
   static LES3DFsd_cState FluxAUSMplus_up_n(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr,
                                            const Vector3D &norm_dir);
   static LES3DFsd_pState Reflect(const LES3DFsd_pState &W,
                                  const Vector3D &norm_dir);
   static LES3DFsd_pState MovingWall(const LES3DFsd_pState &Win,
                                     const LES3DFsd_pState &Wout,
                                     const Vector3D &norm_dir,				 
                                     const Vector3D &wall_velocity,
                                     const Vector3D &pressure_gradient,
                                     const int &TEMPERATURE_BC_FLAG);
   static LES3DFsd_pState NoSlip(const LES3DFsd_pState &Win, 
                                 const LES3DFsd_pState &Wout, 
                                 const Vector3D &norm_dir,  
                                 const Vector3D &pressure_gradient,
                                 const int &TEMPERATURE_BC_FLAG);

   // subfilter stress tensor 
   Tensor3D lambda(const LES3DFsd_pState &dWdx, 
                   const LES3DFsd_pState &dWdy,
                   const LES3DFsd_pState &dWdz,
                   const int &Flow_Type,
                   const double &Volume);
   Tensor3D lambda(const LES3DFsd_pState &dWdx, 
                   const LES3DFsd_pState &dWdy,
                   const LES3DFsd_pState &dWdz,
                   const int &Flow_Type,
                   const double &Volume) const;
   // Heat flux vector 
   Vector3D qflux_t(const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int &Flow_Type,
                    const double &Volume);
   Vector3D qflux_t(const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int &Flow_Type,
                    const double &Volume) const;
   static LES3DFsd_cState Src_t(const LES3DFsd_pState &Wc,
                                const LES3DFsd_pState &dWdx,
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int &Flow_Type,
                                const double &Volume);
   
   /* Index operators */
   double &operator[](int index);
   const double &operator[](int index) const;
   
   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file,
                                const LES3DFsd_pState &W);
   friend istream& operator >> (istream &in_file,
                                LES3DFsd_pState &W);
  
/*  // Destructors */
/*   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; */
/*   specdata = NULL;} */
  
/*   void Deallocate(void){ if(spec != NULL) delete[] spec; */
/*   spec = NULL;  } */
  
/*   ~LES3DFsd_pState(){ Deallocate(); } */
};

/***********************************************************
              LES3D_TheramllyPerfect_cState      
************************************************************/
class LES3DFsd_cState : public NavierStokes3D_ThermallyPerfect_cState {
  public:

   double rhoC;
   double rhoFsd;
   double rhok;
   
// constructors
   
   LES3DFsd_cState(): 
      NavierStokes3D_ThermallyPerfect_cState(){rhoC = rho*ONE; rhoFsd = rho*ONE; rhok = rho*MILLION;}
   
   LES3DFsd_cState(const double &value): 
      NavierStokes3D_ThermallyPerfect_cState(value){rhoC = value; rhoFsd = value; rhok = value;}
   
   LES3DFsd_cState(const double &d, const double &vx,
                   const double &vy, const double &vz,
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En ){rhoC = dC; rhoFsd = dFsd; rhok = dk;}
   
   LES3DFsd_cState(const double &d, const double &vx,
                   const double &vy, const double &vz,
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk,
                   Species *rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){rhoC = dC; rhoFsd = dFsd; rhok = dk;}
   
   LES3DFsd_cState(const double &d, const Vector3D &V,
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk):
      NavierStokes3D_ThermallyPerfect_cState(d, V, En){rhoC = dC; rhoFsd = dFsd; rhok = dk;}
   
   LES3DFsd_cState(const double &d, const Vector3D &V,
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk,
                   const double &frac):
      NavierStokes3D_ThermallyPerfect_cState(d, V, En, frac){rhoC = dC; rhoFsd = dFsd; rhok = dk;}
 
   LES3DFsd_cState(const double &d, const double &vx, 
                   const double &vy, const double &vz, 
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk,
                   const double &rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac){rhoC = dC; rhoFsd = dFsd; rhok = dk;}
   
   LES3DFsd_cState(const double &d, const Vector3D &dV,
                   const double &En, const double &dC, 
                   const double &dFsd, const double &dk,
                   Species *rhomfrac):
      NavierStokes3D_ThermallyPerfect_cState(d, dV, En, rhomfrac){rhoC = dC; rhoFsd = dFsd; rhok = dk;}

//this is needed for the operator overload returns!!!!
   LES3DFsd_cState(const LES3DFsd_cState &U) { rhospec = NULL; rho = DENSITY_STDATM; set_initial_values(); Copy(U); }
                           
   LES3DFsd_pState W(void) ; 
   LES3DFsd_pState W(void) const;
   LES3DFsd_pState W(const LES3DFsd_cState &U) const;

   /* Index operators */
   double &operator[](int index);
   const double &operator[](int index) const;
  
//Read in ns species data
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);
   /* VACUUM **/
   void Vacuum(){ Euler3D_ThermallyPerfect_cState::Vacuum(); rhok = ZERO; rhoC = ZERO; rhoFsd = ZERO;}  

   void Copy(const LES3DFsd_cState &U);
   double T(void) const;
   double a(void) const;     
   static double Mref;               // Mref for Precondtioning (normally set to incoming freestream Mach)
/*    double Laminar_Flame_Speed; */
/*    double Laminar_Flame_Thickness; */
/*    double Adiabatic_Temperature; */
/*    double Equivalence_Ratio; */
/*    double Reactants_Density; */
   double p(void) const;      //pressure
   double C(void) const;
   double Fsd(void) const;
   double k(void) const;
   bool negative_scalarcheck(void);
   bool negative_speccheck(const int &step); //-ve mass frac check and sets small -ve c's to ZERO
   bool Unphysical_Properties_Check(const int n);
  /************** Premixed combustion ************************/
   LES3DFsd_cState premixed_mfrac(const LES3DFsd_pState &Wo);

   /* Binary arithmetic operators. */
   LES3DFsd_cState operator +(const LES3DFsd_cState &U) const;
   LES3DFsd_cState operator -(const LES3DFsd_cState &U) const;
   LES3DFsd_cState operator *(const double &a) const;
   friend LES3DFsd_cState operator *(const double &a, const LES3DFsd_cState &U);
   LES3DFsd_cState operator /(const double &a) const;
   
   double operator *(const LES3DFsd_cState &U) const;
   LES3DFsd_cState operator ^(const LES3DFsd_cState &U) const;
   
  /* Assignment Operator. */
   LES3DFsd_cState& operator =(const LES3DFsd_cState &U);
  /* Shortcut arithmetic operators. */
   LES3DFsd_cState& operator +=(const LES3DFsd_cState &U);
   LES3DFsd_cState& operator -=(const LES3DFsd_cState &U);
      
  /* Unary arithmetic operators. */
   friend LES3DFsd_cState operator -(const LES3DFsd_cState &U);
  
  /* Relational operators. */
   friend int operator ==(const LES3DFsd_cState &U1,
                          const LES3DFsd_cState &U2);
   friend int operator !=(const LES3DFsd_cState &U1,
                          const LES3DFsd_cState &U2);

  /* Input-output operators. */
   friend ostream& operator << (ostream &out_file,
                                const LES3DFsd_cState &U);
   friend istream& operator >> (istream &in_file,
                                LES3DFsd_cState &U);
  
/*    /\* Destructors *\/ */
/*   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata;  */
/*   specdata = NULL;  */
/*   if(Schmidt != NULL) delete[] Schmidt;  */
/*   Schmidt = NULL;  */
/*   } */
/*   void Deallocate(void){ if(rhospec != NULL) delete[] rhospec;  */
/*   rhospec = NULL;  } */
  
/*   ~LES3DFsd_cState(){ Deallocate(); } */
};

inline LES3DFsd_cState LES3DFsd_pState::U(void){
   LES3DFsd_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhoC = rho*C;
   Temp.rhoFsd = rho*Fsd;
   Temp.rhok = rho*k;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
   return Temp;
}

inline LES3DFsd_cState LES3DFsd_pState::U(void)const{
   LES3DFsd_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhoC = rho*C;
   Temp.rhoFsd = rho*Fsd;
   Temp.rhok = rho*k;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
   return Temp;
}

inline LES3DFsd_cState LES3DFsd_pState::U(const LES3DFsd_pState &W){
  if(ns == W.ns){ //check that species are equal
    LES3DFsd_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    Temp.rhoC = W.rho*W.C;
    Temp.rhoFsd = W.rho*W.Fsd;
    Temp.rhok = W.rho*W.k;
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

inline LES3DFsd_cState U(const LES3DFsd_pState &W) {
   LES3DFsd_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   Temp.rhoC = W.rho*W.C;
   Temp.rhoFsd = W.rho*W.Fsd;
   Temp.rhok = W.rho*W.k;
   for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
   }
   return Temp;
}

inline LES3DFsd_pState LES3DFsd_cState::W(void){
   LES3DFsd_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.C = C();
   Temp.Fsd = Fsd();
   Temp.k = k();
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   return Temp;
}

inline void LES3DFsd_pState::Copy(const LES3DFsd_pState &W){
   rho = W.rho;
   v = W.v; 
   p = W.p;  
   C = W.C;
   Fsd = W.Fsd;
   k = W.k;
   for( int i=0; i<ns; i++){
      spec[i] = W.spec[i];
   }
}

//----------------- Index Operator ------------------------/
inline double& LES3DFsd_pState::operator[](int index) {
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
      return C;
   case 7:
      return Fsd;
   case 8:
      return k;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline const double& LES3DFsd_pState::operator[](int index) const {
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
      return C;
   case 7:
      return Fsd;
   case 8:
      return k;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline LES3DFsd_pState LES3DFsd_cState::W(void)const{
   LES3DFsd_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.C = C();  
   Temp.Fsd = Fsd();
   Temp.k = k();  
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   return Temp;
}

inline LES3DFsd_pState LES3DFsd_cState::W(const LES3DFsd_cState &U) const{
  if(ns == U.ns){ //check that species are equal   
    LES3DFsd_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    Temp.p = U.p();
    Temp.C = U.C();  
    Temp.Fsd = U.Fsd();
    Temp.k = U.k();  
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

inline LES3DFsd_pState W(const LES3DFsd_cState &U) {
  LES3DFsd_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  Temp.p = U.p();
  Temp.C = U.C();  
  Temp.Fsd = U.Fsd();
  Temp.k = U.k();
  for(int i=0; i<U.ns; i++){
    Temp.spec[i] = U.rhospec[i]/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  }
  return Temp;
}

inline void LES3DFsd_cState::Copy(const LES3DFsd_cState &U){
   rho = U.rho;
   rhov = U.rhov; 
   E = U.E; 
   rhoC = U.rhoC;
   rhoFsd = U.rhoFsd;
   rhok = U.rhok;
   for( int i=0; i<ns; i++){ 
      rhospec[i] = U.rhospec[i];
   } 
}

//index operators
inline double& LES3DFsd_cState::operator[](int index) {
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
      return rhoC;
   case 7:
      return rhoFsd;
   case 8:
      return rhok;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline const double& LES3DFsd_cState::operator[](int index) const{
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
      return rhoC;
   case 7:
      return rhoFsd;
   case 8:
      return rhok;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline bool LES3DFsd_cState::negative_speccheck(const int &step) {
  double sum(ZERO);
  double temp(ZERO);

  //-------- Negative Check ------------//
  for(int i=0; i<ns-1; ++i){
    temp = rhospec[i].c/rho;
    if(temp > ONE){            //check for > 1.0
      rhospec[i].c =rho;
      temp = ONE;

    } else if(temp < ZERO){   //check for -ve
      if(temp > -1e-8){  //check for small -ve and set to ZERO
	rhospec[i].c = ZERO;
	temp = ZERO;
      } else {

	if( step < 10){
	  return false;
	} else {

	  rhospec[i].c = ZERO;
	  temp = ZERO;
	}
      }
    } else {
      // cout<<"\n negative_speccheck else";
    }
    sum += temp;
  }

//   rhospec[ns-1].c = rho*(ONE - sum); //PUSH error into NS-1 species, probably N2

  temp = max(ONE- sum, ZERO);    //Spread Error across species
  sum += temp;
  rhospec[ns-1].c = rho*temp;
  for(int i=0; i<ns; ++i){
    rhospec[i].c = rhospec[i].c*(ONE/sum);
  }
  
  return true;
}

inline bool LES3DFsd_cState::negative_scalarcheck(void) {
  double LOCAL_TOL = MICRO;
  //-------- Negative Check ------------//
  if ( rhoC < LOCAL_TOL ) { rhoC = ZERO; }
  if ( rhoFsd < LOCAL_TOL ) { rhoFsd = ZERO; }
  if ( rhok < LOCAL_TOL ) { rhok = ZERO; }
  //  if ( rhoC/rho > 0.999 || rhoC/rho < 0.001 ) { rhoFsd = ZERO; }
  return (1);
}

/**************************************************************
  Unphysical_Properties_Check
 ***************************************************************/
inline bool LES3DFsd_cState::Unphysical_Properties_Check(const int n){

  // check for nan's, inf's etc.... debugging !!!
     if (rho <= ZERO || !negative_speccheck(n) || es() <= ZERO || !negative_scalarcheck() ) {
    cout << "\n " << CFFC_Name()
	 << " LESPremixed2D ERROR: Negative Density || Energy || Mass Fractions || Progress Variable || Flame Surface Density : \n"
	 << *this << endl;
    return false;
  }else{
    return true ;
  }
}

#endif //end LES3D_THERMALLYPERFECT_STATE_INCLUDED 

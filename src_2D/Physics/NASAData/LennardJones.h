/************************************************************
 ************************************************************
  File:  LennardJones.h                                    
                                                           
  This header defines the container class the evaluation of 
  transport properties using Lennar-Jones potentials.             
  See the Theory section of the Chemkin manual... i.e.   The 
  trans.dat file was stolen from Cantera.  
                                                           
  Author:  Marc R.J. Charest                               
                                                           
  Date:  Feb 09, 2007                                      
                                                           
 ************************************************************
 ************************************************************/
#ifndef _LENNARD_JONES_INCLUDED
#define _LENNARD_JONES_INCLUDED

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <iomanip>



class LennardJonesData;


/************** MISC. FUNCTIONS *******************************************
  Some necessary functions required to compute the transport properties 
  using Champman-Enskog relations
 *************************************************************************/
// collsion integrals
double Omega11(const double &Ts); 
double Omega22(const double &Ts);
double dOmega22dT(const double &T, const double &eps);

// rotational relaxation function
double F(const double &Ts);      

// Binary Diffusion Coefficient
double BinaryDiff( const LennardJonesData& s1, 
		   const LennardJonesData& s2, 
		   const double &T, const double &P);

//conversion functions (externally defined elsewhere)
extern int stoi(const string);
extern double stof(const string);


/************************************************************
  LennardJonesData Class

  This class contains all the parameters required for the 
  computation of the viscosity, thermal conductivity, and 
  binary diffusion coefficient for a given species.
  See the Theory section of the Chemkin manual for details
  of the computation.
 ************************************************************/

class LennardJonesData {

 public:

  // Lennard-Jones Parameters (only used when tranport_type = TRANSPORT_LENNARD_JONES)
  int g;           // geometry factor ( 0=monatomic, 1=linear, 2=non-linear)
  double eps;      // Lennard-Jones potential well depth [J]
  double sigma;    // Lennard-Jones collision diameter [m]
  double mu;       // dipole moment [m^1.5 J^0.5]
  double alpha;    // polarizability [cubic m]
  double Zrot;     // Rotational relaxation parameter at 298K
  bool   polar;    // true if polar molecule, false if non-polar
  double mol_mass; // molar mass [kg/kmole]

  // default constructor
 LennardJonesData() : g(0), eps(0), sigma(0), mu(0), alpha(0), 
                      Zrot(0), polar(false), mol_mass(0) {}

  // member functions
  void GetData(const char* data_filename, const string spec);
  double Viscosity(const double &T) const;
  double dViscositydT(const double &T) const;
  double ThermalConduct(const double &HeatCapacity_v,
			const double &T) const;
};

/************************************************************
 * Description: Viscosity computed Chapman-Enskog theory    * 
 *              See 'The properties of gases and liquids by *
 *              Reid et al.(1987).  This is Units [N s/m^2] *
 ************************************************************/
inline void LennardJonesData::GetData(const char* datafilename, 
				      const string spec) {
  // declares
  ifstream in;
  string species;
  string line, temp;
  int species_l;
  
  // creat the search string
  species = spec + ' ';
  species_l = species.length();

  //---------------------------------
  // Get transport properties  
  //---------------------------------
 
  // open the transport data file
  in.open(datafilename);
  if(in.fail()){ 
    cerr<<"\nError opening file: tran.dat" <<endl;
    exit(0); 
  }
	   
  // search for the species
  do{
    getline(in,line);
  } 
#ifdef _GNU_GCC_296 //redhat 7.3 gcc-2.96 
  while( line.compare(species,0,species_l) != 0 && !in.eof() );
#else //alpha & redhat 9.0 gcc-3.22
  while( line.compare(0,species_l,species) != 0 && !in.eof() );
#endif

  // if we reached the end of the file and haven't found the species, error
  if(in.eof()){
    cerr<<" Species: "<<species<<" not in tran.dat"<<endl;
    exit(0);
  }

  // get the geometry factor
  g = stoi( line.substr(19,1) );
  
  // Get the potential well-depth / boltzman constant
  // and un-normalize it. [well-depth / boltzman] -> well-depth
  eps = stof( line.substr(20,10) ) * BOLTZMANN;
 
  // get the collision diameter in Angstroms and convert to m
  sigma = stof( line.substr(30,10) ) * 1.0E-10;  
  
  // get the dipole moment in [Debye] and convert to [m^1.5 J^0.5]
  mu = stof( line.substr(40,10) );
  polar = (mu > TOLER); // true if polar (mu>0), false if non-polar (mu=0)
  mu *= 1.0E-25*sqrt(10.0); 
  
  // get the polarizability in [Angstroms^3] and convert to [m^3]
  alpha = stof( line.substr(50,10) ) * 1.0E-30; 
 
  // get the rotational relaxation parameter
  Zrot = stof( line.substr(60,10) );
  
  // close the transport data file
  in.close();

  //cout << "\n Species "<< species <<": g = " << g << "\t eps = " << eps
  //     << "\t sigma = " << sigma << "\n\t mu = " << mu  
  //     << "\t alpha = " << alpha << "\t Zrot = " << Zrot;;

 
}

/************************************************************
 * Description: Viscosity computed Chapman-Enskog theory    * 
 *              See 'The properties of gases and liquids by *
 *              Reid et al.(1987).  This is Units [N s/m^2] *
 ************************************************************/
inline double LennardJonesData::Viscosity(const double &T) const { // [K]
  
  /* Compute */
  double v = (5.0/16.0) * sqrt( PI * mol_mass  * BOLTZMANN * T / AVOGADRO ) 
    / ( PI * sigma*sigma * Omega22( T*BOLTZMANN/eps ) );
  
  return v;
}

inline double LennardJonesData::dViscositydT(const double &T) const { // [K]
  
  /* Compute */
  double prefactor = (5.0/16.0) * sqrt( PI * mol_mass  * BOLTZMANN  / AVOGADRO )
    /( PI * sigma*sigma );
  double v1 = 1.0/( 2.0 * Omega22( T*BOLTZMANN/eps ) * sqrt(T) );
  double v2 = - sqrt(T) * dOmega22dT(T,eps) / ( pow(Omega22(T*BOLTZMANN/eps), 2.0) );

  return prefactor*(v1+v2);
 
}

/************************************************************
 * Description: Thermal conductivity computed Chapman-Enskog*
 *              theory.  See 'The properties of gases and   *
 *              liquids by Reid et al.(1987).               *
 *              Units [W/(m K)]                             *
 *                                                          *
 * Heat Capacity (const Volume):  HeatCapacity_v   J/(kg*K) *
 *                                                          *
 ************************************************************/
inline double LennardJonesData::ThermalConduct(const double &HeatCapacity_v, // [J/(kg*K)]
					       const double &T) const {      // [K]
 
  /* Declares */
  double l;         // thermal conductivity
  double A, B;      // constants
  double Cv_trans;  // translational
  double Cv_rot;    // rotational
  double Cv_vib;    // vibrational
  double f_trans;   // translational
  double f_rot;     // rotational
  double f_vib;     // vibrational
  double Dkk;       // self diffusion coefficient
  double v;         // viscosity
  double d;         // density
  double Z;         // rotational relaxation collision number
  double Tref=298.0;// reference temperature [K]
  double a;
  double eta = Viscosity(T);
  double P = 1.0;   // rho=fn(P) * D=fn(1/P) ==> rho*D independant of P
                    // so we just set it to unity
  
  
  /* Determine molar heat capacity relationships */
  
  Cv_trans = 3.0/2.0;
  
  // for a linear  molecule
  if (g==1) {
    Cv_rot = 1.0;
    Cv_vib = HeatCapacity_v*mol_mass/R_UNIVERSAL - 5.0/2.0;
  }
  // for a non-linear  molecule
  else if (g==2) {
    Cv_rot = 3.0/2.0;
    Cv_vib = HeatCapacity_v*mol_mass/R_UNIVERSAL - 3.0;  
  }

  /* Compute thermal conductivity */
  // for a monatomic gas
  if (g==0) {
    f_trans = 5.0/2.0;
    l = eta * R_UNIVERSAL * f_trans*Cv_trans / mol_mass;
  // for a non-monatonic
  } else {
    Dkk = BinaryDiff(*this,*this,T,P);  // self-diffusion coefficient per unity pressure [D~1/P]
    v = eta;
    d = P*mol_mass / (R_UNIVERSAL*T);   // density per unit pressure [rho~P]
    a = d*Dkk/v;                        // [d*rho !~ P]
    
    Z = Zrot * F(Tref*BOLTZMANN/eps) / F(T*BOLTZMANN/eps);
    A = 5.0/2.0 - a;
    B = Z + (2.0/PI)*( 5.0*Cv_rot/3.0 + a );
    f_trans = (5.0/2.0) * ( 1.0 - (2.0/PI)*(Cv_rot/Cv_trans)*(A/B) );
    f_rot = a * ( 1.0 + (2/PI)*(A/B) );
    f_vib = a;

    l = (v*R_UNIVERSAL/mol_mass) * ( f_trans*Cv_trans + f_rot*Cv_rot + f_vib*Cv_vib  );
  }
  

  return l;
  
}


/************************************************************
 * Function: Omega11                                        *
 *                                                          *
 * Description: Collision integral computed using curvefits *
 *              published by Neufeld (1972).                * 
 *                                                          *
 *              T = temperature                             *
 *              Ts = T*kB/eps                               *
 *              where kB = bolzman constant and             *
 *                    eps = potential well-depth [J]        *
 *                                                          *
 ************************************************************/
inline double Omega11(const double &Ts) {

  // declares
  static int err_cnt(0);

  // constants
  static double A = 1.06036;
  static double B = 0.15610;
  static double C = 0.19300;
  static double D = 0.47635;
  static double E = 1.03587;
  static double F = 1.52996;
  static double G = 1.76474;
  static double H = 3.89411;
  
  // make sure reduced temperature is within range
  if ( (Ts<0.3 || Ts>100.0) && err_cnt<5 ) {
    err_cnt++;
    cerr << "NASARP1311dataclass::Omega11() - Warning: Reduced temperature T*= " 
	 << Ts << " outside applicable range (0.3 < T* < 100)." << endl;
    if (err_cnt==5) cerr << "NASARP1311dataclass::Omega11() - Suppressing error output for now.\n";
    //exit(0);
  }
 
  // compute
  double a = A*pow(Ts,-B) + C*exp(-D*Ts) + E*exp(-F*Ts)
    + G*exp(-H*Ts);
  
  return a;
  
}



/************************************************************
 * Function: Omega22                                        * 
 *                                                          *
 * Description: Collision integral computed using curvefits *
 *              published by Neufeld (1972).                * 
 *                                                          *
 *              T = temperature                             *
 *              Ts = T*kB/eps                               *
 *              where kB = bolzman constant and             *
 *                    eps = potential well-depth [J]        *
 *                                                          *
 ************************************************************/
inline double Omega22(const double &Ts) {

  // declares
  static int err_cnt(0);

  // constants
  static double A =  1.16145;
  static double B =  0.14874;
  static double C =  0.52487;
  static double D =  0.77320;
  static double E =  2.16178;
  static double F =  2.43787;
  static double G =  0.0;
  static double H =  0.0;
  static double R = -6.435E-04;
  static double S =  18.0323;
  static double W = -0.76830;
  static double P =  7.27371;
  
  // make sure reduced temperature is within range
  if ( (Ts<0.3 || Ts>100.0) && err_cnt<5 ) {
    err_cnt++;
    cerr << "NASARP1311dataclass::Omega22() - Warning: Reduced temperature T*= " 
	 << Ts << " outside applicable range (0.3 < T* < 100)." << endl;
    if (err_cnt==5) cerr << "NASARP1311dataclass::Omega22() - Suppressing error output for now.\n";
    //exit(0);
  }
  
  // compute
  double a = A*pow(Ts,-B) + C*exp(-D*Ts) + E*exp(-F*Ts)
    + G*exp(-H*Ts) + R*pow(Ts,B)*sin( S*pow(Ts,W) - P );
  
  return a;
}

 
/************************************************************
 * Function: dOmega22dT                                     *
 *                                                          *
 * Description: Derivative of 2,2 collision integral        *
 *              w.r.t temperature computed using curvefits  *
 *              published by Neufeld (1972).                * 
 *                                                          *
 *              T = temperature                             *
 *              Ts = T*kB/eps                               *
 *              where kB = bolzman constant and             *
 *                    eps = potential well-depth [J]        *
 *                                                          *
 ************************************************************/
inline double dOmega22dT(const double &T, const double &eps) {
  
  // declares
  static int err_cnt(0);

  // constants
  static double A =  1.16145;
  static double B =  0.14874;
  static double C =  0.52487;
  static double D =  0.77320;
  static double E =  2.16178;
  static double F =  2.43787;
  static double G =  0.0;
  static double H =  0.0;
  static double R = -6.435E-04;
  static double S =  18.0323;
  static double W = -0.76830;
  static double P =  7.27371;
  
  double ek = eps/BOLTZMANN;
  double Ts = T/ek;
  
  // make sure reduced temperature is within range
  if ( (Ts<0.3 || Ts>100.0) && err_cnt<5 ) {
    err_cnt++;
    cerr << "Warning: Reduced temperature T*= " << Ts 
	 << " outside applicable range (0.3 < T* < 100)." << endl;
    if (err_cnt==5) cerr << "NASARP1311dataclass::dOmega22dT() - Suppressing error output for now.\n";
    //exit(0);
  }

  // compute
  double a = - A*B*pow(Ts,-B-1.0) - C*D*exp(-D*Ts)
    - E*F*exp(-F*Ts) - G*H*exp(-H*Ts) 
    + R*B*pow(Ts,B-1.0)*sin( S*pow(Ts,W) - P )
    + R*S*W*pow(Ts,B+W-1.0)*cos( S*pow(Ts,W) - P );
 
  //  A*B*pow(Ts,-B-1.0)*pow(Ts,-B)  ->  A*B*pow(Ts,-B-1.0)

  return a/ek;
}


 
/************************************************************
 * Function: F                                              *
 *                                                          *
 * Description: Rotational relaxation function of Parker    *
 *              (1959) and Brau and Jonkman (1970). This    *
 *              function allows for the computation of the  *
 *              rotational relaxation parameter at a        *
 *              arbitrary temperature given the value at the*
 *              reference temperature of 298K.              *
 *              F = Z/Zref                                  *
 *                                                          *
 *              T = temperature                             *
 *              Ts = T*kB/eps                               *
 *              where kB = bolzman constant and             *
 *                    eps = potential well-depth [J]        *
 *                                                          *
 ************************************************************/
inline double F(const double &Ts) {
  
  double f = 1.0 + 0.5*pow(PI,1.5)*pow(Ts,-0.5) 
    + ( 0.25*PI*PI+2.0 )/Ts 
    + pow(PI,1.5)*pow(Ts,-1.5);
  
  return f;
  
}


 
/************************************************************
 * Function: BinaryDiff                                     *
 *                                                          *
 * Description: Diffusivity computed Chapman-Enskog         *
 *              theory.  See 'The properties of gases and   *
 *              liquids by Reid et al.(1987).               *
 *              Units [m^2/s]                               *
 *                                                          *
 ************************************************************/
inline double BinaryDiff( const LennardJonesData& s1, 
			  const LennardJonesData& s2, 
			  const double &T, const double &P) {
  
  /* Declares*/
  const double tol = 1.0E-6;
  double eps;     // effective potential well-depth
  double Ts;      // effective reduced temperature
  double sigma;   // effective collision diameter
  double mu;      // effective dipole moment
  double M;       // effective molar mass
  double d;       // diffusion coefficient
  double alpha_n; // non-polar polarizability
  double mu_p;    // polar dipole moment
  double xi;      // polar correction factor
  double f_eps;   // well depth correction factor
  double f_sigma; // collision diameter correction factor
  const LennardJonesData* sp;    // polar species pointer
  const LennardJonesData* sn;    // non-polar species pointer
  
  
  /* Determine effective properties */
  
  // the well depth
  eps = sqrt(s1.eps*s2.eps);
  
  // the reduced temperature
  Ts = T*BOLTZMANN/eps;
  
  // dipole moment
  mu = sqrt(s1.mu*s2.mu);
  
  // collision diameter
  sigma = 0.5*(s1.sigma + s1.sigma);
  
  
  /* Determine corrections */

  //
  // for both polar or both non-polar
  //
  if ( s1.polar == s2.polar ) {
    //cout << "\nBOTH POLAR OR BOTH NON_POLAR";
    
    f_eps = f_sigma = 1.0;
    
  //
  // for one polar and one non-polar
  //
  } else {
    //cout << "\nPOLAR AND NON_POLAR";

    // determine polar one
    sp = (s1.mu>0.0 ? &s1 : &s2); // polar
    sn = (s1.mu>0.0 ? &s2 : &s1); // nonpolar
    
    // the corrections factors
    alpha_n = sn->alpha/pow(sn->sigma,3.0);
    mu_p = sp->mu/sqrt( sp->eps*pow(sp->sigma,3.0) );
    xi = 1.0 + 0.25*alpha_n*mu_p*sqrt(sp->eps/sn->eps);

    f_eps =  xi*xi;
    f_sigma = pow(xi, -1.0/6.0);

  }
  
  /* apply corrections factors */
  eps *= f_eps;
  sigma *= f_sigma;
  
  
  /* Compute diffusivity */
  
  // effective molar mass
  M = s1.mol_mass*s2.mol_mass / ( AVOGADRO * (s1.mol_mass+s2.mol_mass) );
  
  d = (3.0/16.0) * sqrt( 2.0 * PI * pow(BOLTZMANN * T, 3.0)/M ) 
    / ( PI * P * sigma*sigma * Omega11(Ts) );
  
  return d;
  
}


#endif //_TRANS_INCLUDED

/***********************************************************************
 ***********************************************************************
 **  File: SNBCK.h                                                    **
 **                                                                   **
 **  Description:  This header defines the function prototypes,       **
 **                constants, and any class definitions required for  **
 **                the SNBCK non-grey gas radiation model.            **
 **                                                                   **
 **                The SNBCK model is outlined in Liu et al. (2000).  **
 **                Using the Statistical Narrow Band model parameters **
 **                provided by Soufiani and Taine (1997) and the      **
 **                Correlated-k method of Lacis and Oinas (1991),     **
 **                narrow band spectral intensities may be computed.  **
 **                Essentially, given the temperature, pressure, and  **
 **                species composition (CO, CO2, H2O, soot), the      **
 **                absorbsion coefficient for each band, at each      **
 **                quadrature point, may be computed.  Once the RTE   **
 **                is solved, the spectral intensity accross          **
 **                each band can be averaged using the quadrature.    **
 **                Currently, overlapping bands are treated using the **
 **                various approximations outlined by Liu et          **
 **                al. (2001) (optically thin, optically thick,       **
 **                uncorrelated).  Bands are also                     **
 **                lumped together using the strategy proposed by     **
 **                Liu et al. (1999).  An additional option of        **
 **                precalculating the cummulative distribution        **
 **                function at each wavenumber and quadrature point   **
 **                as a function of temperature is provided.          **
 **                This requires the assumption that the SNB model    **
 **                parameter B is independant of species              **
 **                concentration.  The absorbsion coefficients for    **
 **                each species computed from the curvefits are then  **
 **                treated as uncorrelated (i.e. added together).     **
 **                This procedure is outlined in Liu and Smallwood    **
 **                (2004).                                            **
 **                                                                   **
 ** References:  Soufiani and Taine, "High termperature gas radiative **
 **              property parameters of statistical narrow-band model **
 **              for H2O, CO2 and CO, and correlated-K model for H2O  **
 **              and CO2," International Journal of Heat and Mass     **
 **              Transfer, vol. 40, no. 4, pp. 987-991. 1997.         **
 **                                                                   **
 **              F. Liu, G.J. Smallwood, O.L. Gulder, "Application of **
 **              the statistical narrow-band correlated-k method to   **
 **              low-resolution, spectral intenstity and radiative    **
 **              heat transfer calculations - effects of the          **
 **              quadrature scheme." International Journal of Heat    **
 **              and Mass Transfer, v43, 2000. pp. 3119-3135.         **
 **                                                                   **
 **              F. Liu, G.J. Smallwood, O.L. Gulder, "Application    **
 **              of the statistical narrow-band correlated-k method   **
 **              to non-grey gas radiation in CO2-H2O mixtures:       **
 **              approximate treatment of overlapping bands," in      **
 **              Journal of Quantitative Spectroscopy & Radiative     **
 **              Heat Transfer, v68 (2001). pp 401-417.               **
 **                                                                   **
 **              F. Liu, G.J. Smallwood, O.L. Gulder, "Radiation heat **
 **              transfer calculations using the SNBCK method,"       **
 **              AIAA paper 99-3679 (1999).                           **
 **                                                                   **
 **              F. Liu, G.J. Smallwood, "An efficient approach for   **
 **              implementation of the SNB based correlated-k method  **
 **              and its evaluation."  in                             **
 **              Journal of Quantitative Spectroscopy & Radiative     **
 **              Heat Transfer, v84 (2004). pp 465-475.               **
 **                                                                   **
 **              V. Goutiere, A. Charette, L. Kiss, "Comparative      **
 **              performance of nongray gas modelling techniques,"    **
 **              Numerical Heat Transfer, Part B, 41:3, 361-381. 2002 **
 **                                                                   **
 **  Author:  Marc R.J. Charest                                       **
 **                                                                   **
 **  Date:    Dec 17th, 2006                                          **
 ***********************************************************************
 ***********************************************************************/

/*********************************************************************
 ***************************** INCLUDES ******************************
 *********************************************************************/

#ifndef _SNBCK_INCLUDED
#define _SNBCK_INCLUDED

// Required C++ libraries
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// Required header files
#include "../../Math/Math.h"
#include "../../Math/Vector2D.h"
#include "../../Math/SplineFit.h"
#include "../../Math/Quadrature.h"
#include "../../CFD/CFD.h"
#include "../../MPI/MPI.h"
#include "../Planck.h"


// Use cubic splines or polynomial interpolants
#define SNBCK_USE_SPLINES
//#define SNBCK_USE_POLY

/*********************************************************************
 ***************************** CONSTANTS *****************************
 *********************************************************************/

// unit conversion
#define CM_TO_M         1.0E-2

// file length
#define	INPUT_PARAMETER_LENGTH_RTE2D    128

// flags for quadrature type.  quadrature was taken from : 
// F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass Trans., v43, 2000. pp. 3119-3135.
enum QuadType { GAUSS_LEGENDRE, // Gauss-Legendre quadrature scheme
		GAUSS_LOBATTO };// Gauss-Lobatto quadrature scheme

// flags for treatment at overlapping bands
// proposed by Liu et al. (2001), J Quant Spect & Rad Heat Trans, v68 (2001). pp 401-417
enum BandTreatment { SNBCK_OVERLAP_OPTICALLY_THIN, // optically thin
		     SNBCK_OVERLAP_OPTICALLY_THICK,// optically thick
		     SNBCK_OVERLAP_UNCORRELATED,   // uncorrelated
		     SNBCK_OVERLAP_CORRELATED };   // correlated


// flags for inversion of the cummulative distribution function
// proposed by Liu and Smallwood. (2001), J Quant Spect & Rad Heat Trans, v84 (2004). pp 465-475
enum Inversion { SNBCK_EVAL_ONLINE,    // direct newton solve
		 SNBCK_EVAL_PRECALC }; // precalulated and stored in cubic spline

// Pressure broadening parameter cutoff.  This is the maximum value
// of the pressure broadening term, B, so that overflow is avoided when
// exp(PI*B) is evaluated.
enum { SNBCK_BROADENING_CUTOFF = 200 };

/*********************************************************************
 ********************* FUNCTION PROTOTYPES ***************************
 *********************************************************************/

// return value of cummulative distribution function and its derivative
double g( const double k, const double B, const double S );
double g_prime( const double k, const double B, const double S );
void g_lumped( const double*B, const double*S,
	       const int Nstart, const int Nend, const int Nlump,
	       const double k, double &gg, double &dgg);

// return the value of the absorbsion coefficient for a specified quad point
double AbsorptionCoeffSNBCK( const double g, 
			     const double*B, const double*S,
			     const int Nstart, const int Nend);
double PeakAbsorptionCoeffSNBCK( const double B, const double S );

// plot the cummulative distribution function at a certain band
void OutputDistFunction(  const double B, const double S,
			  const double kmin, const double kmax,
			  const char *file );

// compute gas transmissivity
double TransmissivitySNB( const double B, const double S, const double L );

// compute line of site inensity for an isothermal, homogeneous gas
double LineOfSightIntens( const double I1, const double Ib, 
			  const double k, const double L );
double LineOfSightIntens( const double I1, const double Ib, 
			  const double tau );


// convert fortran "D" notation to "E" notation for scientific notation
double FortranToCppD( string str );

// approximation of error function
double erfc_a(const double x);
double erfc_exp(const double x1, const double x2);


/*********************************************************************
 *********************** CLASS DEFINITIONS ***************************
 *********************************************************************/

class SNBCK;
class EM2C;
class SNBCK_Input_Parameters;


/******************* SNBCK CLASS DEFINTION *********************************
 *                                                                         *
 * Class: EM2C                                                             *
 *                                                                         *
 * Description:  The EM2C class is used to store statistical narrow band   *
 *               parameters of the Soufiani and Taine (1997) dataset,      *
 *               EM2C*.  This code was basically copied from the EM2C      *
 *               fortran program provided in the Appendix of Modest's book,*
 *               "Radiative Heat Transfer," 2nd ed., 2003.                 *
 *                                                                         *
 * Usage:   1 -> Declare object                                            *
 *               >> EM2C SNB;                                              *
 *          2 -> Load the database into memory (given the PATH to the data)*
 *               >> SNB.LoadParams(PATH);                                  *
 *          3a-> Given T, P, and species conc., compute the SNB parameters *
 *               'B' and 'S' for each band, for each gas                   *
 *               >> SNB.ComputeSNB( p, T, xco, xh2o, xco2, xo2 );          *
 *          3b-> Alternatively, to compute gamma at a reference state and  *
 *               generate the SNB model parameter S per unit mole X and    *
 *               per unit atm,                                             *
 *               >> SNB.ComputeRefSNB( p0, T0, xco0, xh2o0, xco20, xo20 ); *
 *          4 -> You can now use S and B as you wish.                      *
 *                                                                         *
 * NOTE: Model parameters are stored based on the following units:         *
 *         length - cm                                                     *
 *         press  - atm                                                    *
 *         temp   - K                                                      *
 ***************************************************************************/
class EM2C{
  
 public:

  //------------------------------------------------
  // Objects
  //------------------------------------------------

  // constants
  static const int NCO     =  48;  // number of CO bands
  static const int NH2O    = 367;  // number of H2O bands
  static const int NCO2    =  96;  // number of CO2 bands
  static const int Npoints =  14;  // number of temperature points
  static const int Nbands  = 367;  // number of bands
  static const double Tmin;        // minimum valid temperature [K]
  static const double Tmax;        // max valid temperature [K]

  // datafile paths and names
  char FileSNBCO [INPUT_PARAMETER_LENGTH_RTE2D]; // CO SNB parameter datafile
  char FileSNBH2O[INPUT_PARAMETER_LENGTH_RTE2D]; // H2O SNB parameter datafile
  char FileSNBCO2[INPUT_PARAMETER_LENGTH_RTE2D]; // CO2 SNB parameter datafile
  char FileSNBWN [INPUT_PARAMETER_LENGTH_RTE2D]; // SNB wavenumber datafile

  // model parameters ( k [cm^-1 atm^-1] and 1/delta [cm])
  double kCO[NCO][Npoints], kH2O[NH2O][Npoints], kCO2[NCO2][Npoints]; // k
  double dCO[NCO][Npoints], dH2O[NH2O][Npoints], dCO2[NCO2][Npoints]; // 1/delta

  // wavenumber and bandwidth
  double WaveNo[Nbands], BandWidth[Nbands];  // [cm^-1]

  // spectral indexes for each species relating narrow band to dataset
  int   iCO[Nbands],  iCO2[Nbands],  iH2O[Nbands];

  // logical indicating whether band is active (=1), or inactive (=0)
  int  liCO[Nbands], liCO2[Nbands], liH2O[Nbands], liMix[Nbands];

  // SNB model parameters
  double B_CO  [Nbands], B_CO2  [Nbands], B_H2O[Nbands], 
         B_Thin[Nbands], B_Thick[Nbands];
  double S_CO  [Nbands], S_CO2  [Nbands], S_H2O[Nbands], 
         S_Thin[Nbands], S_Thick[Nbands];
          

  //------------------------------------------------
  // Member Functions
  //------------------------------------------------

  // load the model constants from the data file
  void LoadParams(const char *PATH);

  // compute SNB parameters
  void ComputeRefSNB( const double p,      // pressure [atm]
		      const double T,      // temperature [K]
		      const double xCO,    // mole fraction oh CO
		      const double xH2O,   // mole fraction oh H2O
		      const double xCO2,   // mole fraction oh CO2
		      const double xO2 );  // mole fraction oh O2
  void ComputeSNB( const double p,      // pressure [atm]
		   const double T,      // temperature [K]
		   const double xCO,    // mole fraction oh CO
		   const double xH2O,   // mole fraction oh H2O
		   const double xCO2,   // mole fraction oh CO2
		   const double xO2 );  // mole fraction oh O2

  // multiply S by arbitrary state
  void MultiplySNB( const double p_0,         // pressure [atm]
		    const double xCO_0,       // mole fraction of CO
		    const double xH2O_0,      // mole fraction of H2O
		    const double xCO2_0 );    // mole fraction of CO2

  // compute planck mean absorbsion coefficient
  void PlanckMean( const double T,   // temperature [K]
		   double &kp_CO,    // planck mean abs. of CO
		   double &kp_H2O,   // planck mean abs. of H2O
		   double &kp_CO2 ); // planck mean abs. of CO2
  
  // compute gas transmissivity
  double Transmissivity( const double L, const int i );

private:
  
  // determine the spectral index for each species
  void FindIndex();
  
  // determine interpolation coefficients
  void Interp( const double T, double &slope, int &index);

};


/******************* SNBCK CLASS DEFINTION *********************************
 *                                                                         *
 * Class: SNBCK                                                            *
 *                                                                         *
 * Description:  The SNBCK class is used to compute narrow band            *
 *               intensities. Given the temperature, pressure,  species    *
 *               composition (CO, CO2, H2O, soot), quadrature rule, the    *
 *               absorbsion coefficient for each band, at each quadrature  *
 *               point, are computed.  The absorbsion coefficients can then*
 *               be used in radiative heat transfer computations.  Band    *
 *               averaged values, say for intensity, are then determined   *
 *               by integrating the parameter in question over the         *
 *               cummulative distribution function.                        *
 *                                                                         *
 * Usage:   1 -> Declare object                                            *
 *               >> SNBCK S;                                               *
 *          2 -> Setup SNBCK class with the necessary parameters           *
 *               a - For online inversion of the cummulative distribution  *
 *                   function:                                             *
 *                   >> S.Setup( quad_type, Nquad, Nlump, Optimize,        *
 *                               CURR_PATH, overlap_model );               *
 *               b - For curvefit representation wrt T of the cummulative  *
 *                   distribution function:                                *
 *                   >> S.Setup( quad_type, Nquad, Nlump, Optimize,        *
 *                               CURR_PATH, p_ref, xco_ref, xh2o_ref,      *
 *                               xco2_ref, xo2_ref, Nint );                *
 *          3 -> For a given state, calculate the absorbsion coefficient   *
 *               for each band, for each quadrature point                  *
 *               >> S.CalculateAbsorb( P, T, xco, xh2o, xco2, xo2, xsoot );*
 *          4 -> Compute the band average of a parameter phi given its     *
 *               values at each quadrature point at band 'v':              *
 *               >> avg = S.BandAverage( phi, v );                         *
 *                                                                         *
 * NOTE: Model parameters are stored based on the following units:         *
 *         length - cm                                                     *
 *         press  - atm                                                    *
 *         temp   - K                                                      *
 ***************************************************************************/
class SNBCK{

  
  //------------------------------------------------
  // Objects
  //------------------------------------------------

public:

  // distribution function inversion evaluation type
  int EvalType;

  // Mixture rule to be used
  int MixType;
  
  //the number of quadrature points at each band
  int Nquad;
  int* nquad;

  // the number of wide, lumped narrow, bands
  int Nbands;

  // the nondimensional wavenumber array (quadrature points)
  double* g;

  // the quadrature weight array
  double *w, **ww;

  // the wavenumber center for each wide band [cm^-1]
  double* WaveNo;
  double* BandWidth;

  // spectral indexes relating wide bands to narrow bands
  int *istart, *iend;

  // EM2C dataset class
  EM2C SNB;

  // flags to indicate if wide band active (>0) or transparent (=0)
  int *iCO, *iCO2, *iH2O, *iMix; 

  // an index relating a 2D array (v,i) to an unfolded 1D array (n)
  int **index;

  // an index relating an unfolded 1D array (n) to 2D array (v,i)
  int *band_index;
  int *quad_index;

private:
  
  // Interpolation parameters required for cubic spline fitting of
  // the absorbsion coefficient wrt temperature.  Uniform spacing in 
  // temperature is used.
  int Ninterp;       // number of tempeature interpolation points
  double Tmin, Tmax; // temperature limits
  double dT;         // stepsize between temperature interpolation points (uniform spacing)
  double *Tn;        // temperature points
  double ***k_CO,    // k wrt temp for interpolation
         ***k_CO2,
         ***k_H2O;
  double ***k2_CO,   // second derivative of k wrt temp for interpolation
         ***k2_CO2,
         ***k2_H2O;

  // temporary storage
  double **k;        // absorbsion coefficient [cm^-1]

  //------------------------------------------------
  // Member Functions
  //------------------------------------------------

public:
  
  // Consutructor 1 - components uninitialized
  SNBCK() : SNB(), Nquad(0), nquad(NULL), Nbands(0),  
            g(NULL), ww(NULL), w(NULL), WaveNo(NULL), 
            BandWidth(NULL), istart(NULL), iend(NULL),
            iCO(NULL), iCO2(NULL), iH2O(NULL), iMix(NULL),
            Tn(NULL), k_CO(NULL), k_CO2(NULL), k_H2O(NULL), 
            k2_CO(NULL), k2_CO2(NULL), k2_H2O(NULL), dT(ZERO),
	    Tmax(ZERO), Tmin(ZERO),
            Ninterp(0), MixType(SNBCK_OVERLAP_OPTICALLY_THIN),
            EvalType(SNBCK_EVAL_ONLINE), index(NULL),
            band_index(NULL), quad_index(NULL), k(NULL) {}
  SNBCK( const SNBCK_Input_Parameters &IP,  // input parameters
	 const char *CFFC_PATH )            // Current path
    : SNB(), Nquad(0), nquad(NULL), Nbands(0),  
      g(NULL), ww(NULL), w(NULL), WaveNo(NULL), 
      BandWidth(NULL), istart(NULL), iend(NULL),
      iCO(NULL), iCO2(NULL), iH2O(NULL), iMix(NULL),
      Tn(NULL), k_CO(NULL), k_CO2(NULL), k_H2O(NULL), 
      k2_CO(NULL), k2_CO2(NULL), k2_H2O(NULL), dT(ZERO),
      Tmax(ZERO), Tmin(ZERO),
      Ninterp(0), MixType(SNBCK_OVERLAP_OPTICALLY_THIN),
      EvalType(SNBCK_EVAL_ONLINE), index(NULL),
      band_index(NULL), quad_index(NULL), k(NULL) 
      { Setup(IP, CFFC_PATH); }
  
  // destructor
  ~SNBCK() {  Deallocate(); }

  // return the total number of variables (bands * quad points)
  int NumVar() { 
    return ( index[Nbands-1][ nquad[Nbands-1]-1 ] + 1 );  
  }

  // return the quadrature weight for a certain point, at a specific band
  double Weight( const int v, const int i ) {
    if (MixType==SNBCK_OVERLAP_CORRELATED) return ww[v][i];
    else return w[i];
  }
  
  // setup the quadrature and load SNB model parameters
  void Setup( const SNBCK_Input_Parameters &IP,  // input parameters
	      const char *CFFC_PATH );           // Current path


  // compute SNB model parameters
  void CalculateAbsorb( const double p,        // pressure [atm]
			const double T,        // temperature [K]
			const double xco,      // mole fraction oh CO
			const double xh2o,     // mole fraction oh H2O
			const double xco2,     // mole fraction oh CO2
			const double o2,       // mole fraction oh O2
			const double xsoot );  // volume fraction oh soot 
  void CalculateAbsorb( const double p,        // pressure [atm]
			const double T,        // temperature [K]
			const double xco,      // mole fraction oh CO
			const double xh2o,     // mole fraction oh H2O
			const double xco2,     // mole fraction oh CO2
			const double o2,       // mole fraction oh O2
			const double xsoot,    // volume fraction oh soot 
			double *ka );          // absorbsion coefficient array  [m^-1]
			          
  // return band averaged property
  double BandAverage( const double *phi, const int v );

  // return total intensity
  double TotalIntensity( const double *phi, const int v );

  // calculate the planck distribution for the gas
  void CalculatePlanck( const double T, double* Ib );

  // calculate the radiation source based on the optically thin approximation
  double RadSourceOptThin( const double p,        // pressure [atm]
			   const double T,        // temperature [K]
			   const double xco,      // mole fraction oh CO
			   const double xh2o,     // mole fraction oh H2O
			   const double xco2,     // mole fraction oh CO2
			   const double o2,       // mole fraction oh O2
			   const double xsoot );  // volume fraction oh soot 

  // calculate the planck mean absorbsion coefficient based on the optically thin approximation
  double PlanckMean( const double p,        // pressure [atm]
		     const double T,        // temperature [K]
		     const double xco,      // mole fraction oh CO
		     const double xh2o,     // mole fraction oh H2O
		     const double xco2,     // mole fraction oh CO2
		     const double o2,       // mole fraction oh O2
		     const double xsoot );  // volume fraction oh soot 

private:

  // setup the quadrature
  void SetupQuad( const int quad_type,     // quadrature type
		  const int quad_points ); // number of quadrature points

  // lump bands together
  void SetupBands( const int Nlump,        // number of narrow bands lumped together
		   const bool optimize );  // attempt to optimize band lumping (Nlump>1)


  // precalculate SNB model parameters at the rerference state - precalculated case
  void PreCalculateAbsorb( const double p_ref,    // pressure [atm]
			   const double xco_ref,  // mole fraction oh CO
			   const double xh2o_ref, // mole fraction oh H2O
			   const double xco2_ref, // mole fraction oh CO2
			   const double xo2_ref,  // mole fraction of O2
			   const int Nint );      // number of interpolation points

  // compute absorbsion coeffient - online iversion
  void CalculateAbsorb_Direct( const double p,        // pressure [atm]
			       const double T,        // temperature [K]
			       const double xco,      // mole fraction oh CO
			       const double xh2o,     // mole fraction oh H2O
			       const double xco2,     // mole fraction oh CO2
			       const double xo2,      // mole fraction oh O2
			       const double fsoot );  // volume fraction oh soot 

  // compute absorbsion coeffient - precalculated interpolation
  void CalculateAbsorb_Interp( const double p,        // pressure [atm]
			       const double T,        // temperature [K]
			       const double xco,      // mole fraction oh CO
			       const double xh2o,     // mole fraction oh H2O
			       const double xco2,     // mole fraction oh CO2
			       const double xo2,      // mole fraction oh O2
			       const double fsoot );  // volume fraction oh soot 

  // allocate and deallocate the arrays
  void AllocateQuad();
  void DeallocateQuad();
  void AllocateBands();
  void DeallocateBands();
  void AllocateInterp();
  void DeallocateInterp();
  void AllocateIndex();
  void DeallocateIndex();
  void AllocateStorage();
  void DeallocateStorage();
  void Deallocate();

};



/******************* SNBCK CLASS DEFINTION *********************************
 *                                                                         *
 * Class: SNBCK_Input_Parameters                                           *
 *                                                                         *
 * Description:  This class is used to store the input parameters for the  *
 *               EM2C and SNBCK classes.  This is used in the Rte2D CFFC   *
 *               module.                                                   *
 *                                                                         *
 ***************************************************************************/
class SNBCK_Input_Parameters{
 private:
 public:

  // Evaluation type (precalculated or online inversion)
  int EvaluationType;

  // Quadrature parameters
  int QuadType;    // Rule.
  int QuadPoints;  // Number of points.
  
  // narrow band lumping parameters
  int LumpedBands;       // Number of narrow bands to lump together.
  bool OptimizedLumping; // Attempt to use the optimized (non-uniform) 
                         // lumping stragegy.

  // treatment type at overalapping bands
  int OverlapModel;
  
  // Precalculated absorbsion coefficients parameters
  int IntPoints;   // number of temperature interpolation points (0 for default)
  double p_ref;    // Reference pressure [atm]
  double xco_ref;  // Reference mole fraction CO
  double xh2o_ref; // Reference mole fraction H2O
  double xco2_ref; // Reference mole fraction CO2
  double xo2_ref;  // Reference mole fraction O2

  //* Default Constructor */
  SNBCK_Input_Parameters() {
    EvaluationType   = SNBCK_EVAL_ONLINE;
    QuadType         = GAUSS_LEGENDRE;
    QuadPoints       = 1;
    LumpedBands      = 36;
    OptimizedLumping = false;
    OverlapModel     = SNBCK_OVERLAP_OPTICALLY_THIN;
    IntPoints        = 0; // use points correspoding to EM2C database (14 points)

    // ref state proposed by Liu and Smallwood (2004).
    p_ref            = 1.0; //[atm]
    xco_ref          = 0.0;
    xh2o_ref         = 0.2;
    xco2_ref         = 0.1;
    xo2_ref          = 0.0;
  };

  // output parameters
  void Output(ostream &out_file) const;

  // input file parser
  int Parse_Next_Input_Control_Parameter(char *code, char *value);

  // MPI communicators
  void Broadcast_Input_Parameters();
#ifdef _MPI_VERSION
  void Broadcast_Input_Parameters(MPI::Intracomm &Communicator,
				  const int Source_Rank);
#endif


  // IO operators  
  friend ostream &operator << (ostream &out_file,
		               const SNBCK_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       SNBCK_Input_Parameters &IP);
};


#endif //_SNBCK_INCLUDED

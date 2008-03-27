/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Scatter.h                                                  **
 **                                                                  **
 ** Description: This file defines some necessary functions for      **
 **              the description and evaluation of the scattering    **
 **              phase function using Legendre polynomials.          **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            12/08/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _SCATTER_INCLUDED
#define _SCATTER_INCLUDED 

/********************************************************
 * Necessary  Constants                                 *
 ********************************************************/
// Scatter Phase Functions
enum Scatter_Models { RTE2D_SCATTER_ISO,  // isotropic scattering
		      RTE2D_SCATTER_F1,   // Kim and Lee (1988)   
		      RTE2D_SCATTER_F2,   // Kim and Lee (1988)
		      RTE2D_SCATTER_F3,   // Kim and Lee (1988) 
		      RTE2D_SCATTER_B1,   // Kim and Lee (1988)
		      RTE2D_SCATTER_B2 }; // Kim and Lee (1988)


/********************************************************
 * Struct needed to integrate the phase function over   *
 * over the solid angle.  Cointains information for     *
 * legendre polynomials.                                *
 ********************************************************/
struct legendre_param {
  double *An; // the expansion coefficient array
  int Mn;     // degree of Legendre polynomial 
};




/********************************************************
 * Routine: Legendre                                    *
 *                                                      *
 * Compute the nth degree legendre polynomial at x      *
 *                                                      *
 ********************************************************/
inline double Legendre( const double &x, const int &n) { 

  // initialize
  double P0(1), P1(x), Pn;

  // 0th degree
  if (n == 0) return P0;

  // 1st degree
  if (n==1) return P1;

  //recursively compute the nth degree
  for (int i=2; i<=n; i++) {
    Pn = ( (2*i-1) * x * P1 - (i-1) * P0 ) / i;
    P0 = P1;
    P1 = Pn;
  }

  // return the value
  return Pn;
  

} 


/********************************************************
 * Routine: PhaseFunc                                   *  
 *                                                      *
 * Setup the constants for the phase function.          *
 * Currently, 2 forward scattering (F2 and F3) and 2    *
 * backward scattering (B2 and B3) phase functions have *
 * been implemented.  See Kim and Lee (1988) for more   *
 * information on these.                                *
 *                                                      *
 ********************************************************/
inline double* PhaseFunc( const int type, int &n) {

  // declares
  double* An;

  //
  // setup the expansion coefficients
  //
  switch (type) {
	    
    //------------------------------------------------
    // for Linear isotropic scattering
    //------------------------------------------------
    case (RTE2D_SCATTER_ISO):
    default:      
      // the degree
      n = 1;
      // create the array and set the constants
      An = new double[n];
      An[0]  = 1.00000;
     break;

    //
    // Forward scattering with the F1 phase function of Kim and Lee (1988)
    case (RTE2D_SCATTER_F1):
      // the degree
      n = 13;
      // create the array and set the constants
      An = new double[n];
      An[0]  = 1.00000;
      An[1]  = 2.53602;
      An[2]  = 3.56549;
      An[3]  = 3.97976;
      An[4]  = 4.00292;
      An[5]  = 3.66401;
      An[6]  = 3.01601;
      An[7]  = 1.23304;
      An[8]  = 1.30351;
      An[9]  = 0.53463;
      An[10] = 0.20136;
      An[11] = 0.05480;
      An[12] = 0.01099;
      break;

    //------------------------------------------------
    // Forward scattering with the F2 phase function of Kim and Lee (1988)
    //------------------------------------------------
    case (RTE2D_SCATTER_F2):     
      // the degree
      n = 9;
      // create the array and set the constants
      An = new double[n];
      An[0] = 1.00000;
      An[1] = 2.00917;
      An[2] = 1.56339;
      An[3] = 0.67407;
      An[4] = 0.22215;
      An[5] = 0.04725;
      An[6] = 0.00671;
      An[7] = 0.00068;
      An[8] = 0.00005;
      break;

    //------------------------------------------------
    // Forward scattering with the F3 phase function of Kim and Lee (1988)
    //------------------------------------------------
    case (RTE2D_SCATTER_F3):
      // the degree
      n = 3;
      // create the array and set the constants
      An = new double[n];
      An[0] = 1.00000;
      An[1] = 1.20000;
      An[2] = 0.50000;
      break;

    //------------------------------------------------
    // Backward scattering with the B1 phase function of Kim and Lee (1988)
    //------------------------------------------------
    case (RTE2D_SCATTER_B1):
      // the degree
      n = 6;
      // create the array and set the constants
      An = new double[n];
      An[0] =  1.00000;
      An[1] = -0.56524;
      An[2] =  0.29783;
      An[3] =  0.08571;
      An[4] =  0.01003;
      An[5] =  0.00063;
      break;

    //------------------------------------------------
    // Backward scattering with the B2 phase function of Kim and Lee (1988)
    //------------------------------------------------
    case (RTE2D_SCATTER_B2):
      // the degree
      n = 3;
      // create the array and set the constants
      An = new double[n];
      An[0] =  1.00000;
      An[1] = -1.20000;
      An[2] =  0.50000;
      break;

  } // endswitch

  // return the array
  return An;

}

/********************************************************
 * Routine: PhaseEval                                   *
 *                                                      *
 * Evaluate the phase function.                         *
 *                                                      *
 ********************************************************/
inline double PhaseEval( const double* An, 
			 const int &n,
			 const double &cos_THETA ) {
  double Phi(ZERO);
  for (int i=0; i<n; i++)
    Phi += An[i]*Legendre( cos_THETA, i );
  return Phi;
}


/********************************************************
 * Routine: NormalizePhase                              *
 *                                                      *
 * Normalizes the phase function so that it's           *
 * integration over the solid angle sums to 4*PI        *
 *                                                      *
 ********************************************************/
inline void NormalizePhase( double ****Phi,    // phase function array (m,l,p,q)
			    const int &Npolar, // number polar angles
			    const int* Nazim,  // number azim angles
			    double **omega )   // control angle element size
{
  // declares
  double g;
  
  //
  // Loop over incoming dirs
  //
  for(int m=0 ; m<Npolar ; m++) 
    for(int l=0 ; l<Nazim[m] ; l++) {
      
      // initialize
      g = 0;
      
      //
      // sum outgoing portion
      //
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  g += Phi[m][l][p][q]*omega[p][q]; 

      // compute the constant
      g /= 4.0*PI;

      // if this is zero, don't divide by zero -> something wrong
      if (fabs(g)<TOLER) {
	cerr << "\nScatter.h::NormalizePhase() - Normalization factor "
	     << "is zero -> something wrong.\n";
	exit(-1);
      } // endif

      // if this is already one, no need to normalize
      if (fabs(g-ONE)<TOLER) continue;

      //
      // normalize
      //
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  Phi[m][l][p][q] = Phi[m][l][p][q]/g;

    }  // end for -in-dirs-

}


#endif //end _SCATTER_INCLUDED 

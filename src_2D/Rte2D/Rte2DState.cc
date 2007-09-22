/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: RTEState.cc                                                **
 **                                                                  **
 ** Description: The radiation state class contains the properties   **
 **              of a gray, absorbing, emitting, anisotropically     **
 **              scattering medium.  Additionally, it contains the   **
 **              state (intensity) of a beam of light passing        **
 **              through the medium.  This file defines the state    ** 
 **              class.                                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#include "Rte2DState.h"


/********************************************************
 * Static member initialization                         *
 ********************************************************/
int         Rte2D_State :: Npolar          = 0;     // number of polar elements
int*        Rte2D_State :: Nazim           = NULL;  // number of azimuthal elements
int         Rte2D_State :: Nband           = 0;     // number of bands (& quad pts)
int***      Rte2D_State :: Index           = NULL;  // Index array
int         Rte2D_State :: NUM_VAR_RTE2D   = 0;     // total number of variables in soln state
double**    Rte2D_State :: mu              = NULL;  // x-dir cosine
double**    Rte2D_State :: eta             = NULL;  // y-dir cosine
double**    Rte2D_State :: xi              = NULL;  // z,azimuthal-dir cosine
double**    Rte2D_State :: omega           = NULL;  // solid angle element size
double*     Rte2D_State :: theta           = NULL;  // polar angle grid points
double**    Rte2D_State :: psi             = NULL;  // azimuthal angle grid points
double*     Rte2D_State :: delta_theta     = NULL;  // polar angle grid size
double**    Rte2D_State :: delta_psi       = NULL;  // azimuthal angle grid size
double***** Rte2D_State :: Phi             = NULL;  // phase function
double      Rte2D_State :: Symmetry_Factor = ONE;   // symmetry multiplyer


/**************************************************************************
 *********************** RTE2D SETUP DIRS FUNCTIONS  **********************
 **************************************************************************/

/********************************************************
 * Function: SetDirsDOM                                 *
 *                                                      *
 * Sets the number of control angles in the polar and   *
 * azimuthal directions.  Also loads the direction      *
 * cosines for use in the DOM formulation of the RTE.   *
 * Two types of quadrature are possible:                *
 *    1. SN quadrature of Carlson and Lathrop,          *
 *       "Discrete ordinates angular quadrature of the  *
 *       neutron tranport equation." NTIS report LA3186.*
 *    2. TN quadrature of Thurgood et al., Transactions *
 *       of the ASME 117 (1995) pp 1068-1070            *
 *                                                      *
 * Note that we are using the FVM formulation of the RTE*
 * so the DOM directions cosines must be multiplied by  *
 * the quadrature weight before hand (ie. mu is really  *
 * w*mu).                                               *
 *                                                      *
 * DOM ONLY WORKS FOR SPACE-MARCHING.                   *
 ********************************************************/
void Rte2D_State :: SetDirsDOM(const int Quad_Type,
			       const int Axisymmetric,
			       const char *CFFC_PATH )
{
  // declares
  string line;
  double fact;
  string quad_str;
  char file[INPUT_PARAMETER_LENGTH_RTE2D]; 

  // get file name
  strcpy(file, CFFC_PATH);
  strcat(file,"/data/SN_QUAD/quad.dat");

  // open the quadrature data file
  ifstream in;
  in.open(file, ios::in);
  if(in.fail()){ 
    cerr<<"\nRte2D_State::SetDirsDOM(): Error opening file: quad.dat" <<endl;
    exit(-1); 
  }	   
    
  //------------------------------------------------
  // Determine the direction cosines for the specified quadrature.
  //------------------------------------------------
  switch (Quad_Type) {
    
    // S2 NONSYMMETRIC from LATHROP and CARLSON
    case RTE2D_DOM_S2:
      quad_str = "DOM_S2";    
      break;

    // S4 SYMMETRIC from LATHROP and CARLSON
    case RTE2D_DOM_S4:
      quad_str = "DOM_S4";    
      break;
	
    // S8 SYMMETRIC from LATHROP and CARLSON
    case RTE2D_DOM_S6:
      quad_str = "DOM_S6";    
      break;

    // S8 SYMMETRIC from LATHROP and CARLSON
    case RTE2D_DOM_S8:
      quad_str = "DOM_S8";    
      break;

    // S12 SYMMETRIC from LATHROP and CARLSON
    case RTE2D_DOM_S12:
      quad_str = "DOM_S12";    
      break;

    // T3 SYMMETRIC from Truelove
    case RTE2D_DOM_T3:
      quad_str = "DOM_T3";    
      break;

    //shouldn't get here
    default:
      cerr << "Rte2D_State::SetDirsDOM(): Invalid value for quadrature type flag." << endl;
      exit(-1);
      break;
      
  //------------------------------------------------
  } //end switch
  //------------------------------------------------


  // search for the quadrature
  do{
    getline(in,line);
  } while( line.compare(0,quad_str.length(),quad_str) != 0 && 
	   !in.eof() );

  // if we reached the end of the file and haven't found the quadrature, error
  if(in.eof()){
    cerr<<"Quadrature: "<<quad_str<<" not in quad.dat"<<endl;
    exit(-1);
  }

  // deallocate first
  DeallocateCosines();

  // Get the number of cosines in the polar direction
  in >> Npolar;

  // get the maximum number of cosines in the azimuthal direction 
  // for each polar direction.  
  AllocateDirs();
  NUM_VAR_RTE2D = 0;
  for (int i=0; i<Npolar; i++) {
    in >> Nazim[i];
    NUM_VAR_RTE2D += Nazim[i];
  }
  NUM_VAR_RTE2D *= Nband;

  // now allocate the cosine arrays
  AllocateCosines(RTE2D_SOLVER_DOM);

  // get the weight multiplication factor
  in >> fact;

  //
  // get the cosines, this is stored in a jagged array
  //
  for (int m=0; m<Npolar; m++)
    for (int l=0; l<Nazim[m]; l++) {      

      // depending upon coordinate type
      switch(Axisymmetric) {
      case AXISYMMETRIC_X:
	in >> mu[m][l];  // radial (x)
	in >> xi[m][l];  // azimuthal
	in >> eta[m][l]; // axial (y)
	break;
      case AXISYMMETRIC_Y:
	in >> eta[m][l]; // radial (y)
	in >> xi[m][l];  // azimuthal
	in >> mu[m][l];  // axial (x)
	break;
      case PLANAR:
      default:
	in >> mu[m][l];  // x
	in >> xi[m][l];  // z
	in >> eta[m][l]; // y
	break;
      } // end switch 

      // the quadrature weight
      in >> omega[m][l];
      omega[m][l] *= fact;

    } // endfor  - dirs
			    
//-------------   DEBUG    ------------------
//   temp = 0;
//   for (int m=0; m<Npolar; m++) 
//     for (int l=0; l<Nazim[m]; l++) 
//       temp += omega[m][l];
//   cout << endl << "Sum Omega/(4 PI) = " << temp/(FOUR*PI) << endl; 
//   cout << endl;
//   for (int m=0; m<Npolar; m++) {
//     for (int l=0; l<Nazim[m]; l++) {
//       cout.width(4);  cout << m+1;
//       cout.width(4);  cout << l+1;
//       cout.width(14); cout << mu[m][l];
//       cout.width(14); cout << eta[m][l];
//       cout.width(14); cout << xi[m][l];
//       cout.width(14); cout << omega[m][l];
//       cout << endl;
//     }
//     cout << endl;
//   }
//-------------   DEBUG    ------------------


  //
  // Premultiply the cosines by the weights 
  // (we are using a FVM formulation).
  //
  for (int m=0; m<Npolar; m++) 
    for (int l=0; l<Nazim[m]; l++) {
      mu[m][l]  *= omega[m][l];
      eta[m][l] *= omega[m][l];
      xi[m][l]  *= omega[m][l];
    } // endfor

  // the symmetry factor (not used, already accounted for in quadrature)
  Symmetry_Factor = ONE;

  // close the file
  in.close();
}

/********************************************************
 * Function: SetDirsFVM                                 *
 *                                                      *
 * Sets the number of control angles in the polar and   *
 * azimuthal directions.  Also computes the direction   *
 * cosines for use in the FVM formulation of the RTE.   *
 * For now, uniform spacing is used.                    *
 ********************************************************/
void Rte2D_State :: SetDirsFVM(const int NumPolarDirs, 
			       const int NumAzimDirs,
			       const int Axisymmetric )
{
  //------------------------------------------------
  // DECLARES
  //------------------------------------------------

  // Plolar angle limits (0 < theta < pi)
  // For 2D cartesion cases, only need theta->[0,pi/2] due to symmetry.
  double theta_min;
  double theta_max;

  // Azimuthal angle limits (0 < psi < 2pi )
  // For 2D axisymmetric cases, only need psi->[0,PI] due to symmetry.
  double psi_min; 
  double psi_max;

  // control element sizes
  double del_psi;
  double del_theta; 
  

  //------------------------------------------------
  // ALLOCATE DIRECTIONS
  //------------------------------------------------

  // deallocate first
  DeallocateCosines();

  // set the total number of dirs in the polar direction
  Npolar = NumPolarDirs;

  // for uniform spacing, set the number of dis in the azimuthal
  AllocateDirs();
  for (int i=0; i<Npolar; i++ ) Nazim[i] = NumAzimDirs;

  // compute the total number of directions x total number of bands
  NUM_VAR_RTE2D = NumAzimDirs * NumPolarDirs * Nband;

  // now allocate the cosine arrays
  AllocateCosines(RTE2D_SOLVER_FVM);
  

  //------------------------------------------------
  // COMPUTE DIRECTION COSINES
  //------------------------------------------------
  // The total solid angle goes from:
  //    0 <= polar angle <= PI,   0 <= azimuthal angle <= 2PI

  // Axisymmetric
  // Note that we go from PI to 2PI (azimuthal) so that the r-direction 
  // cosines are monotonically increasing (Inportant for angular upwinding!).
  // See  Carlson and Lathrop, "Computing Methods in Reactor  Physics", 
  // (1968) pp. 171-266                              *
  if (Axisymmetric) {
    psi_min = PI; 
    psi_max = 2.0*PI;
    theta_min = 0;
    theta_max = PI;
 
  // Cartesian
  // Note that only a half hemisphere has been modelled
  } else {
    psi_min = 0; 
    psi_max = 2.0*PI;
    theta_min = 0;
    theta_max = PI/2.0;
  } // endif

  // the symmetry factor
  Symmetry_Factor = TWO;


  // create the angular grid
  // define the angle elements sizes for uniform spacing
  del_theta = (theta_max-theta_min) / Npolar; 
  theta[0] = theta_min;
  for (int m=0; m<Npolar; m++ ) {

    del_psi = (psi_max-psi_min) / Nazim[m];
    theta[m+1] = theta[m] + del_theta;
    delta_theta[m] = del_theta;

    psi[m][0] = psi_min;
    for (int l=0; l<Nazim[m]; l++ ) {
      psi[m][l+1] = psi[m][l] + del_psi;
      delta_psi[m][l] = del_psi;
    }
  }


  
  // Define the directional vector (a weighted direction cosine).  
  // This vector points in an average direction within the 
  // solid angle element and its length is indicative 
  // of the size of solid angle element size.
  
  //
  // loop over all directions
  //
  for (int m=0; m<Npolar; m++) { // polar direction
    for (int l=0; l<Nazim[m]; l++) { // azimuthal direction
      
      // compute the directional vector components 
      if (Axisymmetric == AXISYMMETRIC_Y) {

	Control_Angle_Avg( theta[m],     // lower polar angle
			   theta[m+1],   // upper polar angle
			   psi[m][l],    // lower azimuthal angle
			   psi[m][l+1],  // upper azimuthal angle
			   eta[m][l],    // radial-direction cosine
			   xi[m][l],     // angular-direction cosine
			   mu[m][l],     // axial-direction cosine
			   omega[m][l] );// solid angle element size

     } else if (Axisymmetric == AXISYMMETRIC_X) {

	Control_Angle_Avg( theta[m],     // lower polar angle
			   theta[m+1],   // upper polar angle
			   psi[m][l],    // lower azimuthal angle
			   psi[m][l+1],  // upper azimuthal angle
			   mu[m][l],     // radial-direction cosine
			   xi[m][l],     // angular-direction cosine
			   eta[m][l],    // axial-direction cosine
			   omega[m][l] );// solid angle element size

      } else {

	Control_Angle_Avg( theta[m],     // lower polar angle
			   theta[m+1],   // upper polar angle
			   psi[m][l],    // lower azimuthal angle
			   psi[m][l+1],  // upper azimuthal angle
			   mu[m][l],     // x-direction cosine
			   eta[m][l],    // y-direction cosine
			   xi[m][l],     // z-direction cosine
			   omega[m][l] );// solid angle element size

      }
      
      // multiply by the symmetry factor
      mu[m][l]    *= Symmetry_Factor;
      eta[m][l]   *= Symmetry_Factor;
      xi[m][l]    *= Symmetry_Factor;
      omega[m][l] *= Symmetry_Factor;
    
    } // endif - azimuthal

  } // endif - polar

//-------------   DEBUG    ------------------
//   temp = 0;
//   for (int m=0; m<Npolar; m++) 
//     for (int l=0; l<Nazim[m]; l++) 
//       temp += omega[m][l];
//   cout << endl << "Sum Omega/(4 PI) = " << temp/(FOUR*PI) << endl; 

//   for (int m=0; m<Npolar; m++) {
//     for (int l=0; l<Nazim[m]; l++) {
//       cout.width(4);  cout << m+1;
//       cout.width(4);  cout << l+1;
//       cout.width(14); cout << mu[m][l];
//       cout.width(14); cout << eta[m][l];
//       cout.width(14); cout << xi[m][l];
//       cout.width(14); cout << omega[m][l];
//       cout << endl;
//     }
//     cout << endl;
//   }
//-------------   DEBUG    ------------------
  
}


/**************************************************************************
 ********************** RTE2D SETUP PHASE FUNCTIONS  **********************
 **************************************************************************/


/********************************************************
 * Function: SetupPhaseDOM                              *
 *                                                      *
 * Sets up the phase function for the DOM.  Isotropic   *
 * and anisptropic scattering may be modelled.  For     *
 * anisotropic scattering, the forward and backward     *
 * scatteing phase functions of Kim and Lee, Int. J.    *
 * Heat MAss Transfer 31 (8) (1988) pp 1711-1721 are    *
 * used.  Note that the phase function is assummed      *
 * independant of wavenumber.                           *
 ********************************************************/
void Rte2D_State :: SetupPhaseDOM( const int type ) {

  //
  // initialize
  //
  double* An;   // the expansion coefficients array
  int Mn;       // the degree
  double g;     // constant for the normalization of the phase function
  Vector3D in_dir, out_dir;  // incoming and outgoing scattered directions
  double dprod;
  int v;

  //------------------------------------------------
  // linear isotropic scattering
  //------------------------------------------------

  // If this is linear isotropic scattering, the phase function
  // is 1 every where.  No need to numerically average it.
  if ( type == RTE2D_SCATTER_ISO ) {

    for (v=0; v<Nband; v++) 
      for(int m=0 ; m<Npolar ; m++) 
	for(int l=0 ; l<Nazim[m] ; l++) 
	  for(int p=0 ; p<Npolar ; p++) 
	    for(int q=0 ; q<Nazim[p] ; q++) 
	      Phi[v][m][l][p][q] = ONE;
    return ;

  } // endif 
  

  //------------------------------------------------
  // anisotropic scatering
  //------------------------------------------------
  // Since the DOM quadrature  will be used to integrate 
  // the phase function over the solid angle,and this will 
  // be performed in the solver, setup is straightforward
  // Note that 3D vectors are used to properly evaluate the 
  // dot product.

  // get the phase function constants
  An = PhaseFunc( type, Mn);

  // For now, we make the phase function independant of wavelength
  v = 0;
  
  //
  // loop over incoming directions
  //
  for(int m=0 ; m<Npolar ; m++) 
    for(int l=0 ; l<Nazim[m] ; l++) {
    
      //
      // Loop over outgoing directions
      //
      for(int p=0 ; p<Npolar ; p++) 
	for(int q=0 ; q<Nazim[p] ; q++) {
	  
	  // incoming ray (divide out omega)
	  in_dir = Vector3D( mu[m][l],eta[m][l],xi[m][l] )/omega[m][l];
	  // outgoing ray (divide out omega)
	  out_dir = Vector3D( mu[p][q],eta[p][q],xi[p][q] )/omega[p][q];
	  // Compute the dot product between the incoming 
	  dprod = dot(in_dir,out_dir);

	  // compute the phase function for the current directions
	  Phi[v][m][l][p][q] = PhaseEval(An, Mn, dprod);
	  
	} //end for out-dirs
    } // end for in-dirs
  
    
  //------------------------------------------------
  // Normalize the phase function.
  //------------------------------------------------
  NormalizePhase( Phi[0],   // phase function array (m,l,p,q)
		  Npolar,   // number polar angles
		  Nazim,    // number azim angles
		  omega );  // control angle element size
    
  //------------------------------------------------
  // now copy the phase function over all bands
  //------------------------------------------------
  for (v=1; v<Nband; v++) 
    for(int m=0 ; m<Npolar ; m++) 
      for(int l=0 ; l<Nazim[m] ; l++) 
	for(int p=0 ; p<Npolar ; p++) 
	  for(int q=0 ; q<Nazim[p] ; q++) 
	    Phi[v][m][l][p][q] = Phi[0][m][l][p][q];


  // clean up memory
  delete[] An;

  return;

}


/********************************************************
 * Function: func                                       *
 *                                                      *
 * This function is needed for integrating the phase    *
 * function.  In SetupPhaseFVM, it is passed to         *
 * SimpsonMultiDim.h:adaptsim() and is integrated.      *
 ********************************************************/
double phase_func(int ndim, double *x, void *params) {

  // declares
  legendre_param P = *(legendre_param *)params;
  double Phi;
  Vector3D in_dir, out_dir;
  double dprod;

  //
  // determine the direction cosines
  //

  // incoming ray
  in_dir.x = sin(x[0])*cos(x[1]);
  in_dir.y = sin(x[0])*sin(x[1]);
  in_dir.z =           cos(x[1]);

  // outgoing ray
  out_dir.x = sin(x[2])*cos(x[3]);
  out_dir.y = sin(x[2])*sin(x[3]);
  out_dir.z =           cos(x[3]);

  // the dot product
  dprod = dot(in_dir,out_dir);

  //
  // compute the phase function for the current directions
  //
  Phi = PhaseEval(P.An, P.Mn, dprod);

  // return the integrand
  return Phi * sin(x[0]) * sin(x[2]);
}


/********************************************************
 * Function: SetupPhaseFVM                              *
 *                                                      *
 * Sets up the phase function for the FVM.  Isotropic   *
 * and anisptropic scattering may be modelled.  For     *
 * anisotropic scattering, the forward and backward     *
 * scatteing phase functions of Kim and Lee, Int. J.    *
 * Heat MAss Transfer 31 (8) (1988) pp 1711-1721 are    *
 * used.  Note that the phase function is assummed      *
 * independant of wavenumber.                           *
 *                                                      *
 * For anisotropic scattering, analytic integration may *
 * be used to determine an average phase function.      *
 * However, this can be computationally intensive       *
 * for complex phase functions.  Analytic integration   *
 * should only be used for accuracy tests. The method   *
 * proposed by Chai, Lee, and Patankar, J. Thermophys   *
 * Heat Transfer 8 (3) (1994) pp. 419-425 is also       *
 * implemented.  It is an approximate evaluation of the *
 * integral using the panel method.                     *
 *                                                      *
 * FIXME - Analytic integration takes too long.  And I  *
 * had some problems obtaining similar results to those *
 * computed using the approximate panel method.  This   *
 * should be checked in more detail.                    *
 *                                                      *
 ********************************************************/
void Rte2D_State :: SetupPhaseFVM( const int type ) {

  //
  // initialize
  //
  // common parameters
  double val;   // value of integration
  double g;     // constant for the normalization of the phase function
  int v;        // band index
  bool analytic_integration(false);  // true->use analytic integration
                                     // false->use discrete approximation

  // panel method integration parameters
  int Nsub(1);  // break each control angle into Nsub x Nsub smaller ones
  double delta_theta_in, delta_theta_out; // incoming / outgoing sub polar angle size
  double delta_psi_in, delta_psi_out;     // incoming / outgoing sub azim angle size
  Vector3D in_dir, out_dir;               // incoming / outgoing direction cosines
  double dprod;                           // dot product
  double omega_in, omega_out;             // solid angle element sizes
  double tsub_lo, tsub_hi;                // polar angle limits
  double psub_lo, psub_hi;                // azimuthal anlge limits

  // Analytic integration parameters
  double* An;   // the polynomial expansion coefficients array
  int Mn;       // the degree of legendre polynomial
  simp_function F;    // function struct for integration
  legendre_param lp;  // function parameters struct for integration
  simp_state S;       // simpson state struct for integration
  simp_params P;      // simpson params struct for integration
  int fevals;         // number of function evaluations
  int err;            // error flag returned by adaptsim

  // setup integration parameters
  P.maxevals = 100000000;
  P.tol = MICRO;

  //------------------------------------------------
  // linear isotropic scattering
  //------------------------------------------------
  // If this is linear isotropic scattering, the phase function
  // is 1 every where.  No need to numerically average it.
  if ( type == RTE2D_SCATTER_ISO ) {

    for (v=0; v<Nband; v++) 
      for(int m=0 ; m<Npolar ; m++) 
	for(int l=0 ; l<Nazim[m] ; l++) 
	  for(int p=0 ; p<Npolar ; p++) 
	    for(int q=0 ; q<Nazim[p] ; q++) 
	      Phi[v][m][l][p][q] = ONE;
    return ;

  } // endif
  

  //================================================
  //
  // anisotropic scatering
  //
  //================================================

  //------------------------------------------------
  // Setup anisotropic scattering by integrating the phase
  // function over the total solid angles.  Use an adaptive simpsons
  // quadrature rule for multidimensional integration
  //------------------------------------------------
  if (analytic_integration) {

    // get the phase function constants and
    // initialize the legendre parameter struct
    An = PhaseFunc( type, Mn);
    lp.An = An;
    lp.Mn = Mn;

    // setup integration ( allocate memory and set parameters )
    malloc_simp_struc( 4, F, S );
    init_simp_struc( F, S );
    
    // setup function struct and params
    F.f = phase_func;
    F.params = &lp;
    
    cout << "\nIntegrating Phase Function Analytically..." << flush;
    
    // For now, we make the phase function independant of wavelength
    v = 0;
    
    //
    // loop over incoming (m,l) and outgoing (p,q) directions
    //
    for(int m=0 ; m<Npolar ; m++) {
      for(int l=0 ; l<Nazim[m] ; l++) {
	
	for(int p=0 ; p<Npolar ; p++) {
	  for(int q=0 ; q<Nazim[p] ; q++) {
	    
	    // initialize before integration
	    init_simp_struc( F, S );
	    
	    // set new integration limits
	    F.xmin[0] = theta[m];  F.xmax[0] = theta[m+1];
	    F.xmin[1] = psi[m][l]; F.xmax[1] = psi[m][l+1];
	    F.xmin[2] = theta[p];  F.xmax[2] = theta[p+1];;
	    F.xmin[3] = psi[p][q]; F.xmax[3] = psi[p][q+1];
	  
	    // compute the integrated phase function
	    err = adaptsim( F, S, P, fevals, val );
	    if (err) { 
	      cerr << "RteState.cc::SetupPhaseFVM() - "
		   << "Error integrating phase function\n";
	      cerr << "Error flag: " << err << endl;
	      exit (-1);
	    }
	    
	    // set
	    // IMPORTANT - We are only modelling a section (typically
	    // half) of the solid anlge range.  Omega hass been weighed 
	    // accordingly.  Thus need to divide out this factor when
	    // averaging phi.
	    Phi[v][m][l][p][q] = val 
 	      / (omega[m][l]/Symmetry_Factor)
 	      / (omega[p][q]/Symmetry_Factor);
	  
	    
	  } //end for -q-
	} // end for -p-
      } // end for -l-
    } // end for -m-
    
    // clean case specific memory
    free_simp_struc( F, S );

  //------------------------------------------------
  // Setup anisotropic scattering by integrating the phase
  // function over the total solid angles.  Use a type
  // of panel method multidimensional integration
  // See Chai et al. (1994).
  //------------------------------------------------
  } else {
    
    // get the phase function constants
    An = PhaseFunc( type, Mn);

    // For now, we make the phase function independant of wavelength
    v = 0;
    
    //
    // loop over incoming (m,l) and outgoing (p,q) directions
    //
    for(int m=0 ; m<Npolar ; m++) {
      for(int l=0 ; l<Nazim[m] ; l++) {
	
	for(int p=0 ; p<Npolar ; p++) {
	  for(int q=0 ; q<Nazim[p] ; q++) {
	    	    
	    // break control angles up into Nsub x Nsub sections
	    delta_theta_in  = ( theta[m+1]  - theta[m]  ) / double(Nsub);
	    delta_psi_in    = ( psi[m][l+1] - psi[m][l] ) / double(Nsub);
	    delta_theta_out = ( theta[p+1]  - theta[p]  ) / double(Nsub);
	    delta_psi_out   = ( psi[p][q+1] - psi[p][q] ) / double(Nsub);

	    //
	    // loop over the sub control angles
	    //
	    val = ZERO;
	    for(int r=0 ; r<Nsub ; r++) // -> theta_sub_in
	      for(int s=0 ; s<Nsub ; s++) // -> psi_sub_in
		for(int t=0 ; t<Nsub ; t++) // -> theta_sub_out
		  for(int u=0 ; u<Nsub ; u++) {// -> psi_sub_out

		    // Compute the control angle element limits
		    // and the incoming and outgoing cosines/element size.
		    // Note: cosines are thos at the element centroid.

		    // incoming
		    tsub_lo = theta[m] + delta_theta_in*r;
		    tsub_hi = theta[m] + delta_theta_in*(r+1);
		    psub_lo = psi[m][l] + delta_psi_in*s;
		    psub_hi = psi[m][l] + delta_psi_in*(s+1);

		    Control_Angle_Ctr( tsub_lo,    // lower polar angle
				       tsub_hi,    // upper polar angle
				       psub_lo,    // lower azimuthal angle
				       psub_hi,    // upper azimuthal angle
				       in_dir.x,   // x-direction cosine
				       in_dir.y,   // y-direction cosine
				       in_dir.z,   // z-direction cosine
				       omega_in ); // solid angle element size

		    // outgoing
		    tsub_lo = theta[p] + delta_theta_out*t;
		    tsub_hi = theta[p] + delta_theta_out*(t+1);
		    psub_lo = psi[p][q] + delta_psi_out*u;
		    psub_hi = psi[p][q] + delta_psi_out*(u+1);

		    Control_Angle_Ctr( tsub_lo,     // lower polar angle
				       tsub_hi,     // upper polar angle
				       psub_lo,     // lower azimuthal angle
				       psub_hi,     // upper azimuthal angle
				       out_dir.x,   // x-direction cosine
				       out_dir.y,   // y-direction cosine
				       out_dir.z,   // z-direction cosine
				       omega_out ); // solid angle element size

		    // the dot product
		    dprod = dot(in_dir,out_dir);

		    // compute the phase function
		    val += PhaseEval(An, Mn, dprod)*omega_in*omega_out;
		    
		  } // endfor - Nsub

	    
	    // set
	    // IMPORTANT - We are only modelling a section (typically
	    // half) of the solid anlge range.  Omega hass been weighed 
	    // accordingly.  Thus need to divide out this factor when
	    // averaging phi.
	    Phi[v][m][l][p][q] = val 
	      / (omega[m][l]/Symmetry_Factor)
	      / (omega[p][q]/Symmetry_Factor);
	  
	    
	  } //end for -q-
	} // end for -p-
      } // end for -l-
    } // end for -m-

  //================================================
  } //endif - integration type
  //================================================


  //------------------------------------------------
  // Normalize the phase function.
  //------------------------------------------------
  NormalizePhase( Phi[0],   // phase function array (m,l,p,q)
		  Npolar,   // number polar angles
		  Nazim,    // number azim angles
		  omega );  // control angle element size
  
  //------------------------------------------------
  // now copy the phase function over all bands
  //------------------------------------------------
  for (v=1; v<Nband; v++) 
    for(int m=0 ; m<Npolar ; m++) 
      for(int l=0 ; l<Nazim[m] ; l++) 
	for(int p=0 ; p<Npolar ; p++) 
	  for(int q=0 ; q<Nazim[p] ; q++) 
	    Phi[v][m][l][p][q] = Phi[0][m][l][p][q];


  // clean up memory
  delete[] An;

  return;


}

/**************************************************************************
 ************************* RTE2D SETUP FUNCTIONS  *************************
 **************************************************************************/


/********************************************************
 * Rte2D_State::SetupART_DOM                            *
 *                                                      *
 * To evaluating the term                               *
 *     -(eta/r) * dI/dpsi                               *
 * using finite differences is not strickly             *
 * conservative.  It would result in a strong unphysical*
 * coupling between the intensities in neighboring.     *
 * directions that causes difficulties obtaining a      *
 * converged solution.  Instead, the treatment          *
 * of Carlson and Lathrop (1968) is used.  They used    *
 * neutron conservation to derive a angular             *
 * redistribution coeff. This function evaluates this   *
 * term for each of the directions.  Here we            *
 * pre-evaluate the terms and they are stored           *
 * temporarily during computations in ths state class.  *
 * We have adopted the formulation of Kim and Baek      *
 * (1998) which is applicable for unstructured and      *
 * body fitted meshes.                                  *
 * See:                                                 *
 *   Carlson and Lathrop, "Computing Methods in         *
 *   Reactor Physics", (1968) pp. 171-266               *
 * and                                                  *
 *   Kim and Baek, J Thermophys Heat Transfer 12        *
 *   (1998)  pp 596-599                                 *
 ********************************************************/
void Rte2D_State :: SetupART_DOM( const Vector2D &nfaceE, const double &AfaceE,
				  const Vector2D &nfaceW, const double &AfaceW,
				  const Vector2D &nfaceN, const double &AfaceN,
				  const Vector2D &nfaceS, const double &AfaceS ) {

  // allocate the ART array
  if (alpha==NULL) {
    alpha = new double*[Npolar];
    for (int i=0; i<Npolar; i++) alpha[i] = new double[ Nazim[i]+1 ];	
  } /* endif */

  // Allocate storage for intensity in the special directions.
  // These directions are those with no angular distribution.
  if (I_half==NULL) {
    I_half = new double*[Npolar];
    for (int i=0; i<Npolar; i++) I_half[i] = new double[ Nazim[i]+1 ];	
  } /* endif */

  //
  // loop over directions
  //
  for (int m=0; m<Npolar; m++) {
    alpha[m][0] = ZERO; // ART=0 for first azimuthal direction
    for (int l=0; l<Nazim[m]; l++) {

      alpha[m][l+1] = alpha[m][l] +
	     (nfaceE.x*mu[m][l] + nfaceE.y*eta[m][l])*AfaceE +
	     (nfaceW.x*mu[m][l] + nfaceW.y*eta[m][l])*AfaceW +
	     (nfaceN.x*mu[m][l] + nfaceN.y*eta[m][l])*AfaceN +
	     (nfaceS.x*mu[m][l] + nfaceS.y*eta[m][l])*AfaceS;
      
    } // endfor 
  } // endfor


}




/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************** EXTERNAL FUNCTIONS **************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/


 /**************************************************************************
  *************************** RIEMANN FUNCTIONS  ***************************
  **************************************************************************/


/**********************************************************
 * Routine: Riemann (Riemann solver, n-direction)         *
 *                                                        *
 * This function returns a solution to Riemann problem    *
 * for the 2D advection-diffusion equations in the        *
 * x-direction, returning the intermediate state          *
 * variables along the ray x/t=0. Since the RTE is a      *
 * linear hyperbolic equation, and the intensity in the   *
 * individual directions                                  *
 * are only coupled through the source terms, this is     *
 * basically just upwinding.                              *
 **********************************************************/
Rte2D_State Riemann_n(const Rte2D_State &Ul,
		      const Rte2D_State &Ur,
		      const Vector2D &norm_dir ) {
       
  // declares
  double cos_angle, sin_angle, cosine;
  Rte2D_State Um;

  // determine direction cosines for the frame rotation
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  //
  // Compute the Rieman state
  //
  for (int v=0; v<Ul.Nband; v++) 
    for (int m=0; m<Ul.Npolar; m++) 
      for (int l=0; l<Ul.Nazim[m]; l++) {

	// compute the direction cosine
	cosine = Ul.mu[m][l]*cos_angle + Ul.eta[m][l]*sin_angle;
	  
	// upwind
	if (cosine >= ZERO) {
	  Um.In(v,m,l) = Ul.In(v,m,l);
	} else {
	  Um.In(v,m,l) = Ur.In(v,m,l);
	} // endif 

      } // endfor
  
  return (Um);

}




/*********************************************************
 * Routine: Flux_n (solution flux function, n-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.                       *
 *                                                       *
 *********************************************************/
Rte2D_State Flux_n(const Rte2D_State &Ul,
		   const Rte2D_State &Ur,
		   const Vector2D &norm_dir) {

  Rte2D_State Um, Flux;
  Um = Riemann_n(Ul, Ur, norm_dir);
  Flux = Um.Fn(norm_dir);
  return (Flux);

}

 /**************************************************************************
  **************************** BC FUNCTIONS  *******************************
  **************************************************************************/


/*********************************************************
 * Routine: Gray_Wall                                    *
 *                                                       *
 * This function returns the gray wall solution state in *
 * a given direction given the solution variables and    *
 * the unit normal vector in the direction of interest.  *
 * This is for the time marching version of the code.    *
 *                                                       *
 * The wall is assumed diffuse, ie. the intensity leaving*
 * the wall is the same in all directions, and gray,     *
 * ie. the spectral emittance is the same for all        *
 * wavelengths.  The for our case the outgoing intensity *
 * I_w = \epsilon_w*Ib_w +                               *
 *       (1-\epsilon_w)/PI *                             *
 *       \int_{in} I_w(s)|\vec{n} \dot \vec{s}| {d\Omega}*
 *                                                       *
 * See Chapter 9 of "Radiative Heat Transfer" by Modest  *
 * (2003) for a definition of this term.                 *
 *    G = \int_{4\pi} I(r,s) {d\Omega}                   *
 *                                                       *
 *                                                       *
 * Note that control overhang                            *
 * can occur for surfaces not aligned with the principle *
 * directions.  There are two approaches:                *
 *   A - Do nothing.  ie. if the mean direction of a     *
 *       specified control angle is reflected into       *
 *       another control angle in which ovrehang occurs, *
 *       assume it is wholyw reflected into that control *
 *       angle.                                          *
 *   B - The overhanging control angle is                *
 *       pixelated, ie. the control angle in question is *
 *       subdivided into a number of pixels.  Here the   *
 *       incoming and outgoing portions of the           *
 *       overhanging control angle are treated           *
 *       differently. This is more accurate.             *  
 * Currently we only implelent approach A.               *
 *                                                       *
 * See Murthy and Mathur, Num Heat Transfer, Part B, 33  *
 * (1998) 397-416.                                       * 
 *                                                       *
 *                                                       *
 * TODO: Implement pixelation                            *
 *********************************************************/
Rte2D_State Gray_Wall(const Rte2D_State &U, 
		      const Vector2D &norm_dir, 
		      const double &wall_temperature,
		      const double &wall_emissivity ) {

  double In, dot, temp, sum, wn;
  Rte2D_State Uwall(U);
  
  //
  // loop over each band
  //
  for (int v=0; v<Uwall.Nband; v++) {
    
    //------------------------------------------------
    // For a black wall, the blackbody intensity
    //------------------------------------------------
    if (Medium2D_State::Absorb_Type == MEDIUM2D_ABSORB_GRAY) {
      In = wall_emissivity * BlackBody(wall_temperature);
    } else if (Medium2D_State::Absorb_Type == MEDIUM2D_ABSORB_SNBCK) {
      wn = Medium2D_State::SNBCKdata->WaveNo[ Medium2D_State::SNBCKdata->band_index[v] ];
      In = wall_emissivity * BlackBody(wall_temperature, wn);
      // Note: band_index[v] relates 1D Rte2D_State(v) array to 2D SNBCK(v,i) array
    } else {
      cerr << "\nRte2D_State.cc::Gray_Wall(): Invalid value for absorbsion type flag.\n";
      exit(-1);
    } // endif


    //------------------------------------------------
    // For grey wall.
    //------------------------------------------------
    if ( fabs(1.0-wall_emissivity)>MICRO ) {
      
      // Iterate, summing the first moment over the half range for all
      // outgoing directions.  Note that for consistency, sum should
      // add to pi.  This is only possible for boundaries oriented in
      // principle directions.  For arbitrary boundaries, see Modest
      // (2003).  This formulation is used mainly to try and maintain
      // consistency, escpecially for the FVM.  For Outgoing directions,
      // nx*sx<0 or ny*sy<0
      sum = ZERO;
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  sum = sum + max( ZERO, -dot );
	} // endfor
      
      // Iterate about each direction and add the absorbed radiation.
      // Add this components energy if this direction is towards the
      // wall (i.e. nx*sx>0 or ny*sy>0)
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  temp = ((1.0-wall_emissivity)/sum) * Uwall.In(v,m,l) * dot; 
	  In +=  max( ZERO, temp  );
	} // end for 
      
    //------------------------------------------------
    } // endif
    //------------------------------------------------
    
    //
    // Set the intensity of all outgoing directions.
    // Iterate through each direction.
    //
    for (int m=0; m<Uwall.Npolar; m++)
      for (int l=0; l<Uwall.Nazim[m]; l++) {
	
	dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	
	// We are only concerned with directions that are leaving the wall and 
	// entering the domain.  (nx*sx<0 or ny*sy<0)
	if ( dot < ZERO ) Uwall.In(v,m,l) = In;  
	
      }// endfor
    
    
  } // endfor - bands

  return (Uwall);

}

/*********************************************************
 * Routine: Gray_Wall_Space_March                        *
 *                                                       *
 * This function sets the wall boundary state for the    *
 * space-marching case.                                  *
 *                                                       *
 * The wall is assumed diffuse, ie. the intensity leaving*
 * the wall is the same in all directions, and gray,     *
 * ie. the spectral emittance is the same for all        *
 * wavelengths.  The for our case the outgoing intensity *
 * I_w = \epsilon_w*Ib_w +                               *
 *       (1-\epsilon_w)/PI *                             *
 *       \int_{in} I_w(s)|\vec{n} \dot \vec{s}| {d\Omega}*
 *                                                       *
 * See Chapter 9 of "Radiative Heat Transfer" by Modest  *
 * (2003) for a definition of this term.                 *
 *    G = \int_{4\pi} I(r,s) {d\Omega}                   *
 *                                                       *
 *                                                       *
 * Note that control overhang                            *
 * can occur for surfaces not aligned with the principle *
 * directions.  There are two approaches:                *
 *   A - Do nothing.  ie. if the mean direction of a     *
 *       specified control angle is reflected into       *
 *       another control angle in which ovrehang occurs, *
 *       assume it is wholyw reflected into that control *
 *       angle.                                          *
 *   B - The overhanging control angle is                *
 *       pixelated, ie. the control angle in question is *
 *       subdivided into a number of pixels.  Here the   *
 *       incoming and outgoing portions of the           *
 *       overhanging control angle are treated           *
 *       differently. This is more accurate.             *  
 * Currently we only implelent approach A.               *
 *                                                       *
 * See Murthy and Mathur, Num Heat Transfer, Part B, 33  *
 * (1998) 397-416.                                       * 
 *                                                      *
 *                                                       *
 * TODO: Implement pixelation                           *
 *********************************************************/
void Gray_Wall_Space_March(Rte2D_State &Uwall, 
			   const Vector2D &norm_dir, 
			   const double &wall_temperature,
			   const double &wall_emissivity ) {
  // declares
  double In, dot, temp, sum, wn;

  //
  // loop over each band
  //
  for (int v=0; v<Uwall.Nband; v++) {
    
    //------------------------------------------------
    // for a black wall
    //------------------------------------------------
    if (Medium2D_State::Absorb_Type == MEDIUM2D_ABSORB_GRAY)
      In = wall_emissivity * BlackBody(wall_temperature);
    else if (Medium2D_State::Absorb_Type == MEDIUM2D_ABSORB_SNBCK) {
      wn = Medium2D_State::SNBCKdata->WaveNo[ Medium2D_State::SNBCKdata->band_index[v] ];
      In = wall_emissivity * BlackBody(wall_temperature, wn);
      // Note: band_index[v] relates 1D Rte2D_State(v) array to 2D SNBCK(v,i) array
    } else {
      cerr << "\nRte2D_State.cc::Gray_Wall_Space_March(): Invalid "
	   << "value for absorbsion type flag.\n";
      exit(-1);
    } // endif


    //------------------------------------------------
    // For grey wall.
    //------------------------------------------------
    if ( fabs(1.0-wall_emissivity)>MICRO ) {
      
      // Iterate, summing the first moment over the half range for all
      // outgoing directions.  Note that for consistency, sum should
      // add to pi.  This is only possible for boundaries oriented in
      // principle directions.  For arbitrary boundaries, see Modest
      // (2003).  This formulation is used mainly to try and maintain
      // consistency, escpecially for the FVM.  For Outgoing directions,
      // nx*sx<0 or ny*sy<0
      sum = ZERO;
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  sum = sum + max( ZERO, -dot );
	} // endfor 
      
      // Iterate about each direction and add the absorbed radiation.
      // Add this components energy if this direction is towards the
      // wall (i.e. nx*sx>0 or ny*sy>0)
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  temp = ((1.0-wall_emissivity)/sum) * Uwall.In(v,m,l) * dot;
	  In += max( ZERO, temp  );
	} // end for
      
    //------------------------------------------------
    } // endif 
    //------------------------------------------------
    
    
    //
    // Set the intensity of all outgoing directions.
    // Iterate through each direction.
    //
    for (int m=0; m<Uwall.Npolar; m++)
      for (int l=0; l<Uwall.Nazim[m]; l++) {
	
	dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	
	// We are only concerned with directions that are leaving the wall and 
	// entering the domain.  
	if ( dot < ZERO ) Uwall.In(v,m,l) = In;  
	
      } // endfor
    
    
  } // endfor Nband


}




/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 * The reflective surfaces are assumed smooth and thus  *
 * reflection is specular. Note that control overhang   *
 * can occur for surfaces not aligned with the principle*
 * directions.  There are two approaches:               *
 *   A - Do nothing.  ie. if the mean direction of a    *
 *       specified control angle is reflected into      *
 *       another control angle in which ovrehang occurs,*
 *       assume it is wholyw reflected into that control*
 *       angle.                                         *
 *   B - The overhanging control angle is               *
 *       pixelated, ie. the control angle in question is*
 *       subdivided into a number of pixels.  Here the  *
 *       incoming and outgoing portions of the          *
 *       overhanging control angle are treated          *
 *       differently. This is more accurate.            *  
 * Currently we only implelent approach A.              *
 *                                                      *
 * See Murthy and Mathur, Num Heat Transfer, Part B, 33 *
 * (1998) 397-416.                                      *
 *                                                      *
 *                                                      *
 * TODO: Implement pixelation                           *
 ********************************************************/
Rte2D_State Reflect(const Rte2D_State &U, const Vector2D &norm_dir) {

  //declares
  Vector3D out_dir, dir, in_dir;
  double cos_angle, sin_angle;
  double dcn, dct;
  Rte2D_State Temp(U);  Temp.Zero();
  bool exact_match;
  double cos_phi, dotp;
  int mm, ll;
	    

  //------------------------------------------------
  // Apply the frame rotation and calculate the primitive
  // solution state variables in the local rotated frame
  // defined by the unit normal vector. */
  //------------------------------------------------
  for (int m=0; m<U.Npolar; m++) 
    for (int l=0; l<U.Nazim[m]; l++) {

      // set the incoming direction
      in_dir.x = U.mu[m][l];
      in_dir.y = U.eta[m][l];
      in_dir.z = U.xi[m][l];
	
      // Determine the direction cosine's for the frame
      // rotation.
      dcn =   in_dir.x*norm_dir.x + in_dir.y*norm_dir.y; // normal
      dct = - in_dir.x*norm_dir.y + in_dir.y*norm_dir.x; // tangent
	  
      // Reflect the normal component in the rotated frame.
      dcn = - dcn;
	  
      // Rotate back to the original Cartesian reference frame.
      out_dir.x = dcn*norm_dir.x - dct*norm_dir.y;
      out_dir.y = dcn*norm_dir.y + dct*norm_dir.x;
      out_dir.z = in_dir.z;

      //
      // Now search for the closest OUTGOING direction,
      // i.e. the closest one to out_dir
      //
      exact_match = false;
      dotp = 0; mm = 0; ll = 0;
      for (int p=0; p<U.Npolar && !exact_match; p++) 
	for (int q=0; q<U.Nazim[p] && !exact_match; q++) {

	  // the direction to check
	  dir.x = U.mu[p][q];
	  dir.y = U.eta[p][q];	  	      
	  dir.z = U.xi[p][q];	  	      
	      
	  // compute the angle between the two vectors
	  cos_phi = fabs( dot(dir,out_dir) / (abs(dir)*abs(out_dir)) );

	  // we found an exact match
	  if ( fabs(cos_phi-ONE) < TOLER ) {
	    mm = p;
	    ll = q;
	    exact_match = true;

	  // else, check to see if this is one of the closest
	  // directions.  If it is, add it to the list.
	  } else if (fabs(cos_phi-ONE)<fabs(dotp-ONE)) {
	    mm = p;
	    ll = q;
	    dotp = cos_phi;
	  } // endif - exact match

	} // endfor - outgoing dirs
	  
      //
      // Allert user to inexact match
      //
      if (!exact_match)
	cout << endl
	     << "Warning - Rte2DState.cc::Reflect - Control angle overlap "
	     << "detected.  Add pixelation at the boundaries."
	     << endl
	     << "m = " << setw(4) << m << ", l = " << setw(4) << l << endl
	     << "m' = " << setw(4) << mm << ", l' = " << setw(4) << ll << endl
	     << "mu = " << setw(18) << U.mu[m][l] << ", eta = " << setw(18) << U.eta[m][l] << endl
	     << "mu' = " << setw(18) << U.mu[mm][ll] << ", eta' = " << setw(18) << U.eta[mm][ll];
	  
      //
      // set the reflected intensity
      // Remember: reflection directions invarient with wavenumber
      //
      for (int v=0; v<U.Nband; v++) 
	Temp.In(v,mm,ll) = U.In(v,m,l);

      
  //------------------------------------------------
    } // endfor - incoming dirs
  //------------------------------------------------

  
  /* Return the reflected state. */
  return (Temp);
}


/********************************************************
 * Routine: Reflect_Space_March                         *
 *                                                      *
 * This function sets the "wall" boundary state for the *
 * space-marching case.                                 *
 *                                                      *
 * The reflective surfaces are assumed smooth and thus  *
 * reflection is specular. Note that control overhang   *
 * can occur for surfaces not aligned with the principle*
 * directions.  There are two approaches:               *
 *   A - Do nothing.  ie. if the mean direction of a    *
 *       specified control angle is reflected into      *
 *       another control angle in which ovrehang occurs,*
 *       assume it is wholyw reflected into that control*
 *       angle.                                         *
 *   B - The overhanging control angle is               *
 *       pixelated, ie. the control angle in question is*
 *       subdivided into a number of pixels.  Here the  *
 *       incoming and outgoing portions of the          *
 *       overhanging control angle are treated          *
 *       differently. This is more accurate.            *  
 * Currently we only implelent approach A.              *
 *                                                      *
 * See Murthy and Mathur, Num Heat Transfer, Part B, 33 *
 * (1998) 397-416.                                      *
 *                                                      *
 *                                                      *
 * TODO: Implement pixelation                           *
 ********************************************************/
void Reflect_Space_March(Rte2D_State &U, const Vector2D &norm_dir) {

  Vector3D out_dir, dir;
  Vector2D in_dir;
  double cos_angle, sin_angle;
  double dcn, dct, dot_prod;
  bool exact_match;
  double cos_phi, dotp;
  int mm, ll;
	    

  //------------------------------------------------
  // Apply the frame rotation and calculate the primitive
  // solution state variables in the local rotated frame
  // defined by the unit normal vector. */
  //------------------------------------------------
  for (int m=0; m<U.Npolar; m++) 
    for (int l=0; l<U.Nazim[m]; l++) {

      // set the incoming direction
      in_dir.x = U.mu[m][l];
      in_dir.y = U.eta[m][l];
	
      // don't reflect any directions leaving the reflection plane
      // and entering the domain
      if ( dot(in_dir,norm_dir) < ZERO ) continue;
	  
      // Determine the direction cosine's for the frame
      // rotation.
      dcn =   in_dir.x*norm_dir.x + in_dir.y*norm_dir.y; // normal
      dct = - in_dir.x*norm_dir.y + in_dir.y*norm_dir.x; // tangent
	  
      /* Reflect the normal component in the rotated frame. */
      dcn = - dcn;
	  
      /* Rotate back to the original Cartesian reference frame. */
      out_dir.x = dcn*norm_dir.x - dct*norm_dir.y;
      out_dir.y = dcn*norm_dir.y + dct*norm_dir.x;
      out_dir.z = U.xi[m][l];


      //
      // Now search for the closest OUTGOING direction,
      // i.e. the closest one to out_dir
      //
      exact_match = false;
      dotp = 0;
      for (int p=0; p<U.Npolar && !exact_match; p++) 
	for (int q=0; q<U.Nazim[p] && !exact_match; q++) {

	  // compute the dot product between the direction to check
	  // and the normal direction to ascertain if it is ougoing or
	  // incoming
	  dot_prod = norm_dir.x * U.mu[p][q] + norm_dir.y * U.eta[p][q];
	      
	  // don't look at any directions leaving the domain and
	  // incident on the reflection plane
	  if ( dot_prod > ZERO ) continue;
	      
	  // the direction to check
	  dir.x = U.mu[p][q];
	  dir.y = U.eta[p][q];	  	      
	  dir.z = U.xi[p][q];	  	      
	      
	  // compute the angle between the two vectors
	  cos_phi = fabs( dot(dir,out_dir) / (abs(dir)*abs(out_dir)) );

	  // we found an exact match
	  if ( fabs(cos_phi-ONE) < TOLER ) {
	    mm = p;
	    ll = q;
	    exact_match = true;

	  // else, check to see if this is one of the closest
	  // directions.  If it is, add it to the list.
	  } else if (fabs(cos_phi-ONE)<fabs(dotp-ONE)) {
	    mm = p;
	    ll = q;
	    dotp = cos_phi;
	  } // endif - exact match

	} // endfor - outgoing dirs
	  
      //
      // Allert user to inexact match
      //
      if (!exact_match)
	cout << endl
	     << "Warning - Rte2DState.cc::Reflect_Space_March - Control angle overlap "
	     << "detected.  Add pixelation at the boundaries.";
	  
      //
      // set the reflected intensity
      // Remember: directions invarient with wavenumber
      //
      for (int v=0; v<U.Nband; v++) 
	U.In(v,mm,ll) = U.In(v,m,l);

      
  //------------------------------------------------
    } // endfor - incoming dirs
  //------------------------------------------------

}


/**************************************************************************
 ******************** RESTRICT/PROLONG FUNCTIONS  *************************
 **************************************************************************/


/********************************************************
 * Routine: Restrict_NonSol                             *
 *                                                      *
 * This function restricts coefficients that are not    *
 * part of the solution state from a fine grid to a     *
 * coarse one.  This is required for the Multigrid      *
 * specializations.                                     *
 *                                                      *
 ********************************************************/
/********************************************************
COMMENTED OUT FOR NOW

void Restrict_NonSol( Rte2D_State &Uc, const double &Ac,          // coarse cell state, volume
		      const Rte2D_State &Uf1, const double &Af1,  // fine cell state, volume
		      const Rte2D_State &Uf2, const double &Af2,  // fine cell state, volume
		      const Rte2D_State &Uf3, const double &Af3,  // fine cell state, volume
		      const Rte2D_State &Uf4, const double &Af4 ) // fine cell state, volume
{
  //
  // Loop over each band
  //
  for ( int v=0; v<Uc.Nband; v++) {
    
    Uc.kappa[v] = 
      ( Uf1.kappa[v] * Af1 + Uf2.kappa[v] * Af2 +
        Uf3.kappa[v] * Af3 + Uf4.kappa[v] * Af4 ) / Ac;	
	  
    Uc.sigma[v] = 
      ( Uf1.sigma[v] * Af1 + Uf2.sigma[v] * Af2 +
        Uf3.sigma[v] * Af3 + Uf4.sigma[v] * Af4 ) / Ac;		  
    
    Uc.Ib[v] = ( Uf1.Ib[v] * Af1 + Uf2.Ib[v] * Af2 +
		 Uf3.Ib[v] * Af3 + Uf4.Ib[v] * Af4 ) / Ac;		  

  } // endfor Nband 

  
}

********************************************************/


/********************************************************
 * Routine: Restrict_NonSol _Boundary_Ref_States        *
 *                                                      *
 * This function restricts coefficients that are not    *
 * part of the solution state from a fine grid to a     *
 * coarse one.  This is required for the Multigrid      *
 * specializations.                                     *
 *                                                      *
 ********************************************************/
/********************************************************
COMMENTED OUT FOR NOW

void Restrict_NonSol_Boundary_Ref_States( 
		      Rte2D_State &Uc,                              // coarse cell state, volume
		      const Rte2D_State &Uf_l, const double &Af_l,  // left fine cell state, volume
		      const Rte2D_State &Uf_r, const double &Af_r ) // right fine cell state, volume
{
  //
  // Loop over each band
  //
  for ( int v=0; v<Uc.Nband; v++) {
    
    Uc.kappa[v] = 
      ( Uf_l.kappa[v] * Af_l + Uf_r.kappa[v] * Af_r ) / ( Af_r+Af_l );	
	  
    Uc.sigma[v] = 
      ( Uf_l.sigma[v] * Af_l + Uf_r.sigma[v] * Af_r ) / ( Af_r+Af_l );	
    
    Uc.Ib[v] = ( Uf_l.Ib[v] * Af_l + Uf_r.Ib[v] * Af_r ) / ( Af_r+Af_l );
	
  } // endfor Nband 


}

********************************************************/


/********************************************************
 * Routine: Prolong_NonSol                              *
 *                                                      *
 * This function prolongs coefficients that are not     *
 * part of the solution state from a fine grid to a     *
 * coarse one.  This is required for the Multigrid      *
 * specializations. Bilinear interpolation is used.     *
 *                                                      *
 ********************************************************/
/********************************************************
COMMENTED OUT FOR NOW

int Prolong_NonSol( const Rte2D_State &Uc1, const Vector2D &Xc1,  // fine cell state, cell center
		    const Rte2D_State &Uc2, const Vector2D &Xc2,  // fine cell state, cell center
		    const Rte2D_State &Uc3, const Vector2D &Xc3,  // fine cell state, cell center
		    const Rte2D_State &Uc4, const Vector2D &Xc4,  // fine cell state, cell center
		    const Vector2D XcP, Rte2D_State &Uf )         // coarse cell state, cell center
{
  int tmp;
  int error_flag = 0;

  //
  // Loop over each band
  //
  for ( int v=0; v<Uf.Nband && error_flag==0; v++) {
    
    error_flag = Bilinear_Interpolation( Uc1.kappa[v], Xc1, Uc2.kappa[v], Xc2, 
					 Uc3.kappa[v], Xc3, Uc4.kappa[v], Xc4,
					 XcP, Uf.kappa[v] );

    if (error_flag == 0)
      error_flag = Bilinear_Interpolation( Uc1.sigma[v], Xc1, Uc2.sigma[v], Xc2, 
					   Uc3.sigma[v], Xc3, Uc4.sigma[v], Xc4,
					   XcP, Uf.sigma[v] );
    
    if (error_flag == 0)
      error_flag = Bilinear_Interpolation( Uc1.Ib[v], Xc1, Uc2.Ib[v], Xc2, 
					   Uc3.Ib[v], Xc3, Uc4.Ib[v], Xc4,
					   XcP, Uf.Ib[v] );
  } // endfor Nband 


  return error_flag;
}


********************************************************/


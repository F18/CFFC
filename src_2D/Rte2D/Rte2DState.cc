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

#ifndef _RTE2D_STATE_INCLUDED
#include "Rte2DState.h"
#endif // _RTE2D_STATE_INCLUDED   


/********************************************************
 * Static member initialization                         *
 ********************************************************/
int         Rte2D_State :: Npolar          = 0;
int*        Rte2D_State :: Nazim           = NULL;
int         Rte2D_State :: Nband           = 0;
int         Rte2D_State :: Ntot            = 0;
int***      Rte2D_State :: Index           = NULL;
int         Rte2D_State :: NUM_VAR_RTE2D   = 0;
double**    Rte2D_State :: mu              = NULL;
double**    Rte2D_State :: eta             = NULL;
double**    Rte2D_State :: xi              = NULL;
double**    Rte2D_State :: omega           = NULL;
double*     Rte2D_State :: theta           = NULL;
double**    Rte2D_State :: psi             = NULL;
double*     Rte2D_State :: delta_theta     = NULL;
double**    Rte2D_State :: delta_psi       = NULL;
double***** Rte2D_State :: Phi             = NULL;
int         Rte2D_State :: Symmetry_Factor = 1;
int         Rte2D_State :: RTE_Type        = RTE2D_SOLVER_FVM;


/**************************************************************************
 *********************** RTE2D SETUP DIRS FUNCTIONS  **********************
 **************************************************************************/

/********************************************************
 * Function: SetDirs                                    *
 *                                                      *
 * Wrapper function for setting the number of control   *
 * angles in the polar and                              *
 * azimuthal directions.  Also loads the direction      *
 * cosines for use in the DOM formulation of the RTE.   *
 ********************************************************/
void Rte2D_State :: SetDirs(const int NumPolarDirs, 
			    const int NumAzimDirs,
			    const int Quad_Type,
			    const int Axisymmetric )
{
  switch(RTE_Type) {
  case RTE2D_SOLVER_FVM:
    SetDirsFVM(NumPolarDirs, NumAzimDirs, Quad_Type, Axisymmetric );
    break;
  case RTE2D_SOLVER_DOM:
    SetDirsDOM(NumPolarDirs, NumAzimDirs, Quad_Type, Axisymmetric );
    break;
  default:
    cerr << "Rte2D_State::SetDirs - Invalid flag for RTE solver\n";
    exit(1);
  } /* endswitch */
  
}

/********************************************************
 * Function: SetDirsDOM                                 *
 *                                                      *
 * Sets the number of control angles in the polar and   *
 * azimuthal directions.  Also loads the direction      *
 * cosines for use in the DOM formulation of the RTE.   *
 * DOM only works for space marching.                   *
 ********************************************************/
void Rte2D_State :: SetDirsDOM(const int NumPolarDirs, 
			       const int NumAzimDirs,
			       const int Quad_Type,
			       const int Axisymmetric )
{
  // declares
  string line;
  double fact;
  string quad_str;

  // open the quadrature data file
  ifstream in;
  in.open("Rte2D/quad.dat", ios::in);
  if(in.fail()){ 
    cerr<<"\nRte2D_State::SetDirsDOM(): Error opening file: quad.dat" <<endl;
    exit(-1); 
  }	   
    
  // Determine the direction cosines for the specified quadrature.
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
      
  } //end switch


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
  Ntot = 0;
  for (int i=0; i<Npolar; i++) {
    in >> Nazim[i];
    Ntot += Nazim[i];
  }

  // set the total number of directions x total number of bands
  NUM_VAR_RTE2D = Ntot;
  
  // now allocate the cosine arrays
  AllocateCosines();

  // get the weight multiplication factor
  in >> fact;

  // get the cosines, this is stored in a jagged array
  for (int m=0; m<Npolar; m++) {
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
      default:
	in >> mu[m][l];  // x
	in >> xi[m][l];  // z
	in >> eta[m][l]; // y
	break;
      }/* end switch */

      in >> omega[m][l];
      omega[m][l] *= fact;


    } /* endfor */
  } /* endfor */
			    
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


  // Premultiply the cosines by the weights 
  // (we are using a FVM formulation).
  for (int m=0; m<Npolar; m++) 
    for (int l=0; l<Nazim[m]; l++) {
      mu[m][l]  *= omega[m][l];
      eta[m][l] *= omega[m][l];
      xi[m][l]  *= omega[m][l];
    } /* endfor */

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
			       const int Quad_Type,
			       const int Axisymmetric )
{

  /* DECLARES */

  // a temporary storage variable
  double temp, temp1, temp2;
  
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
  
 
  /* ALLOCATE DIRECTIONS */

  // deallocate first
  DeallocateCosines();

  // set the total number of dirs in the polar direction
  Npolar = NumPolarDirs;

  // for uniform spacing, set the number of dis in the azimuthal
  AllocateDirs();
  for (int i=0; i<Npolar; i++ ) Nazim[i] = NumAzimDirs;

  // compute the total number of directions x total number of bands
  Ntot = NumAzimDirs * NumPolarDirs * Nband;
  NUM_VAR_RTE2D = Ntot;
  
  // now allocate the cosine arrays
  AllocateCosines();
  


  /* COMPUTE DIRECTION COSINES */

  // Axisymmetric
  // Note that we go from PI to 2PI (azimuthal) so that the r-direction 
  // cosines are monotonically increasing (Inportant for angular upwinding!).
  if (Axisymmetric) {
    psi_min = PI; 
    psi_max = 2*PI;
    theta_min = 0;
    theta_max = PI;
 
  // Cartesian
  // Note that only a half hemisphere has been modelled
  } else {
    psi_min = 0; 
    psi_max = 2*PI;
    theta_min = 0;
    theta_max = PI/2;
  }

  // the symmetry factor
  Symmetry_Factor = 2;


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
  
  // polar direction
  for (int m=0; m<Npolar; m++) {
    
    
    // azimuthal direction
    for (int l=0; l<Nazim[m]; l++) { 
      
      
      // integrate over the polar angle
      temp1 = 0.5*(theta[m+1]-theta[m]) 
	- 0.25*( sin(2*theta[m+1]) - sin(2*theta[m]) );
      temp2 = -0.25*( cos(2.0*theta[m+1]) - cos(2.0*theta[m]) );

      // compute the directional vector components 
      if (Axisymmetric == AXISYMMETRIC_Y) {
	mu[m][l]  = temp2 * ( psi[m][l+1] - psi[m][l]  );             //axial
	eta[m][l]  = temp1 * ( sin(psi[m][l+1]) - sin(psi[m][l])  );  //radial
 	xi[m][l] = - temp1 * ( cos(psi[m][l+1]) - cos(psi[m][l])  );  //angular
     } else if (Axisymmetric == AXISYMMETRIC_X) {
	mu[m][l]  = temp1 * ( sin(psi[m][l+1]) - sin(psi[m][l])  );   //radial
	eta[m][l]  = temp2 * ( psi[m][l+1] - psi[m][l]  );            //axial
	xi[m][l] = - temp1 * ( cos(psi[m][l+1]) - cos(psi[m][l])  );  //angular
      } else {
	mu[m][l]  = temp1 * ( sin(psi[m][l+1]) - sin(psi[m][l])  );   //x
	eta[m][l] = - temp1 * ( cos(psi[m][l+1]) - cos(psi[m][l])  ); //y
	xi[m][l]  = temp2 * ( psi[m][l+1] - psi[m][l]  );             //z
      }
      
      // solid angles
      omega[m][l] = - ( cos(theta[m+1]) - cos(theta[m]) ) 
	* ( psi[m][l+1] - psi[m][l] );

      // multiply by the symmetry factor
      mu[m][l]    *= Symmetry_Factor;
      eta[m][l]   *= Symmetry_Factor;
      xi[m][l]    *= Symmetry_Factor;
      omega[m][l] *= Symmetry_Factor;
    
    }

  }

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
 * Function: SetupPhase                                 *
 *                                                      *
 * Wrapper function for setting up the phase function   *
 * for use in the  either formulation of the RTE.       *
 ********************************************************/
void Rte2D_State :: SetupPhase( const int type )
{
  switch(RTE_Type) {
  case RTE2D_SOLVER_FVM:
    SetupPhaseFVM( type );
    break;
  case RTE2D_SOLVER_DOM:
    SetupPhaseDOM( type );
    break;
  default:
    cerr << "Rte2D_State::SetupPhase - Invalid flag for RTE solver\n";
    exit(1);
  } /* endswitch */
  
}



/********************************************************
 * Function: SetupPhaseDOM                              *
 *                                                      *
 * Sets up the phase function for                       *
 * use in the DOM formulation of the RTE.               *
 ********************************************************/
void Rte2D_State :: SetupPhaseDOM( const int type ) {

  // initialize
  double* An;   // the expansion coefficients array
  int Mn;       // the degree
  double g;     // constant for the normalization of the phase
		// function
  Vector3D in_dir, out_dir;  // incoming and outgoing scattered directions
  double dprod;
  int v;

  // if this is linear isotropic scattering, the phase function
  // is 1 every where.  No need to numerically average it.
  if ( type == RTE2D_SCATTER_ISO ) {

    for (v=0; v<Nband; v++) 
      for(int m=0 ; m<Npolar ; m++) 
	for(int l=0 ; l<Nazim[m] ; l++) 
	  for(int p=0 ; p<Npolar ; p++) 
	    for(int q=0 ; q<Nazim[p] ; q++) 
	      Phi[v][m][l][p][q] = ONE;
    return ;

  } /* endif */
  

  // else we are dealing with anisotropic scatering

  // get the phase function constants
  An = PhaseFunc( type, Mn);


  // For now, we make the phase function independant of wavelength
  v = 0;
  
  // Setup linear anisotropic scattering.  Since the DOM quadrature
  // will be used to integrate the phase function over the solid angle,
  // and this will be performed in the solver, setup is straightforward
  for(int m=0 ; m<Npolar ; m++) {
    for(int l=0 ; l<Nazim[m] ; l++) {
      
      for(int p=0 ; p<Npolar ; p++) {
	for(int q=0 ; q<Nazim[p] ; q++) {
	  
	  // initialize
	  Phi[v][m][l][p][q] = ZERO;

	  // incoming ray
	  in_dir = Vector3D( mu[m][l],eta[m][l],xi[m][l] );
	  // outgoing ray
	  out_dir = Vector3D( mu[p][q],eta[p][q],xi[p][q] );
	  // Compute the dot product between the incoming 
	  dprod = dot(in_dir,out_dir) / abs(in_dir) / abs(out_dir);

	  // compute the phase function for the current directions
	  for (int i=0; i<Mn; i++)
	    Phi[v][m][l][p][q] += An[i]*Legendre( dprod, i );

	  
	} //end for -q-
      } // end for -p-
    } // end for -l-
  } // end for -m-
  
    
  // Normalize the phase function.
  // First, compute the factor
  for(int m=0 ; m<Npolar ; m++) {
    for(int l=0 ; l<Nazim[m] ; l++) {
      
      // initialize
      g = 0;
      
      // sum
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  g += Phi[v][m][l][p][q]*omega[p][q]; 
      
      // if this is zero, don't divide by zero
      if (g<TOLER) continue;
      
      // compute the constant
      g = g/(4.0*PI);

      // if this is already one, no need to normalize
      if (fabs(g-ONE)<TOLER) continue;

      // normalize
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  Phi[v][m][l][p][q] = Phi[v][m][l][p][q]/g;
	
    }  // end for -l-
  } // end for -m-
  
  // now copy the phase function over all bands
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
 * function.  It is passed to                           *
 * SimpsonMultiDim.h:adaptsim()                         *
 ********************************************************/
double phase_func(int ndim, double *x, void *params) {

  // declares
  legendre_param P = *(legendre_param *)params;
  double Phi = 0;
  Vector3D in_dir, out_dir;  // incoming and outgoing scattered directions
  double dprod;

  // determine the direction cosines

  // incoming ray
  in_dir.x = sin(x[0])*cos(x[1]);
  in_dir.y = sin(x[0])*sin(x[1]);
  in_dir.z =           cos(x[1]);

  // outgoing ray
  out_dir.x = sin(x[2])*cos(x[3]);
  out_dir.y = sin(x[2])*sin(x[3]);
  out_dir.z =           cos(x[3]);

  // the dot product
  dprod = dot(in_dir,out_dir) / abs(in_dir) / abs(out_dir);

  // compute 
  for (int i=0; i<P.Mn; i++)
    Phi += P.An[i]*Legendre( dprod, i );

  // return the integrand
  return Phi * sin(x[0]) * sin(x[2]);
}


/********************************************************
 * Function: SetupPhaseFVM                              *
 *                                                      *
 * Sets up the solid angle averaged phase function for  *
 * use in the FVM formulation of the RTE.               *
 ********************************************************/
void Rte2D_State :: SetupPhaseFVM( const int type ) {

  // initialize
  double* An;   // the expansion coefficients array
  int Mn;       // the degree
  double g;     // constant for the normalization of the phase
		// function
  simp_function F;    // function struct for integration
  legendre_param lp;  // function parameters struct for integration
  simp_state S;       // state struct for integration
  simp_params P;      // params struct for integration
  int fevals;         // number of function evaluations
  double val;         // value of integration
  int err;            // error flag returned by adaptsim
  int v;


  // if this is linear isotropic scattering, the phase function
  // is 1 every where.  No need to numerically average it.
  if ( type == RTE2D_SCATTER_ISO ) {

    for (v=0; v<Nband; v++) 
      for(int m=0 ; m<Npolar ; m++) 
	for(int l=0 ; l<Nazim[m] ; l++) 
	  for(int p=0 ; p<Npolar ; p++) 
	    for(int q=0 ; q<Nazim[p] ; q++) 
	      Phi[v][m][l][p][q] = ONE;
    return ;

  } /* endif */
  

  // else we are dealing with anisotropic scatering

  // get the phase function constants
  An = PhaseFunc( type, Mn);
  lp.An = An;
  lp.Mn = Mn;

  // setup integration ( allocate memory and set parameters )
  malloc_simp_struc( 4, F, S );
  init_simp_struc( F, S );
  
  // setup function
  F.f = phase_func;
  F.params = &lp;

  // setup integration parameters
  P.maxevals = 100000000;
  P.tol = 1e-6;

  cout << "\nIntegrating Phase Funcion..." << flush;
  

  // For now, we make the phase function independant of wavelength
  v = 0;
  
  // Setup linear anisotropic scattering by integrating the phase
  // function over the total solid angles.  Use an adaptive simpsons
  // quadrature rule for multidimensional integration
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
	    cerr << "RteState.cc::SetupPhase - "
		 << "Error integrating phase function\n";
	    cerr << "Error flag: " << err << endl;
	    exit (-1);
	  }
	  
	  // set
	  Phi[v][m][l][p][q] = val 
	    / (omega[m][l]/double(Symmetry_Factor))
	    / (omega[p][q]/double(Symmetry_Factor));
	  // IMPORTANT - THE FACTOR 4 COMES FROM SYMMETRY
	  // WE ARE ONLY MODELLING HALF THE SOLID ANGLE RANGE,
	  // THUS OMEGA HAS BEEN MULTIPLIED BY 2.  NEED TO
	  // DIVIDE OUT THIS FACTOR WHEN AVERAGING PHI.
	  
	  
	} //end for -q-
      } // end for -p-
    } // end for -l-
  } // end for -m-
  
    
  // Normalize the phase function.
  // First, compute the factor
  for(int m=0 ; m<Npolar ; m++) {
    for(int l=0 ; l<Nazim[m] ; l++) {
      
      // initialize
      g = 0;
      
	// sum
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  g += Phi[v][m][l][p][q]*omega[p][q]; 
      
      // if this is zero, don't divide by zero
      if (g<TOLER) continue;
      
      // compute the constant
      g = g/(4.0*PI);

      // if this is already one, no need to normalize
      if (fabs(g-ONE)<TOLER) continue;

      // normalize
      for(int p=0 ; p<Npolar ; p++)
	for(int q=0 ; q<Nazim[p] ; q++) 
	  Phi[v][m][l][p][q] = Phi[v][m][l][p][q]/g;
	
    }  // end for -l-
  } // end for -m-
  
  // now copy the phase function over all bands
  for (v=1; v<Nband; v++) 
    for(int m=0 ; m<Npolar ; m++) 
      for(int l=0 ; l<Nazim[m] ; l++) 
	for(int p=0 ; p<Npolar ; p++) 
	  for(int q=0 ; q<Nazim[p] ; q++) 
	    Phi[v][m][l][p][q] = Phi[0][m][l][p][q];

  cout << " norm = " << g << "...done.\n" << flush;


  // clean up memory
  delete[] An;
  free_simp_struc( F, S );

  return;


}

/**************************************************************************
 ************************* RTE2D SETUP FUNCTIONS  *************************
 **************************************************************************/

/********************************************************
 * Function: SetGas                                     *
 *                                                      *
 * Sets the number of gas bands to discretize the entire*
 * spectrum into.                                       *
 ********************************************************/
void Rte2D_State :: SetGas(  )
{
  Nband = 1;
}



/********************************************************
 * Rte2D_State::SetupART_DOM                            *
 *                                                      *
 * Setup axisymmetric angular redistribution            *
 * coefficients as proposed by Carlson and Lathrop      *
 * (1968).  This is for the DOM space march version.    *
 ********************************************************/
void Rte2D_State :: SetupART_DOM( const Vector2D nfaceE, const double AfaceE,
				  const Vector2D nfaceW, const double AfaceW,
				  const Vector2D nfaceN, const double AfaceN,
				  const Vector2D nfaceS, const double AfaceS ) {

  double temp;

  // allocate the ART array
  if (alpha==NULL) {
    alpha = new double*[Npolar];
    for (int i=0; i<Npolar; i++) alpha[i] = new double[ Nazim[i]+1 ];	
  } /* endif */
  if (I_half==NULL) {
    I_half = new double*[Npolar];
    for (int i=0; i<Npolar; i++) I_half[i] = new double[ Nazim[i]+1 ];	
  } /* endif */

  // loop over directions
  for (int m=0; m<Npolar; m++) {
    alpha[m][0] = ZERO; // ART=0 for first azimuthal direction
    for (int l=0; l<Nazim[m]; l++) {

      temp = (nfaceE.x*mu[m][l] + nfaceE.y*eta[m][l])*AfaceE +
	     (nfaceW.x*mu[m][l] + nfaceW.y*eta[m][l])*AfaceW +
	     (nfaceN.x*mu[m][l] + nfaceN.y*eta[m][l])*AfaceN +
	     (nfaceS.x*mu[m][l] + nfaceS.y*eta[m][l])*AfaceS;
      alpha[m][l+1] = alpha[m][l] + temp;
      
    } /* endfor */
  } /* endfor */


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
 * variables along the ray x/t=0.                         *
 *                                                        *
 **********************************************************/
Rte2D_State Riemann_n(const Rte2D_State &Ul,
		      const Rte2D_State &Ur,
		      const Vector2D &norm_dir ) {
       
    double cos_angle, sin_angle, cosine;
    Rte2D_State Um;

    /* determine direction cosines for the frame rotation */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;


    /* Copy the state properties */

    for (int v=0; v<Ul.Nband; v++) {
      Um.Ib[v] = HALF*(Ul.Ib[v]+Ur.Ib[v]);
      Um.kappa[v] = HALF*(Ul.kappa[v]+Ur.kappa[v]);
      Um.sigma[v] = HALF*(Ul.sigma[v]+Ur.sigma[v]);
    }


    /* Compute the Rieman state */
    
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
	  } /* endif */

	} /* endfor*/

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
  Flux.zero_non_sol();
  
  return (Flux);


}

 /**************************************************************************
  **************************** BC FUNCTIONS  *******************************
  **************************************************************************/


/*********************************************************
 * Routine: Gray_Wall                                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.                       *
 *                                                       *
 *********************************************************/
Rte2D_State Gray_Wall(const Rte2D_State &U, 
		      const Vector2D &norm_dir, 
		      const double &wall_temperature,
		      const double &wall_emissivity ) {

  double In, dot, temp, sum;
  Rte2D_State Uwall(U);
  
 
  // The intensity leaving the wall is the same for all outgoing
  // intensities

  // loop over each band
  for (int v=0; v<Uwall.Nband; v++) {
    
    // for a black wall
    In = wall_emissivity * Ib(wall_temperature);

    // For grey wall.
    if ( fabs(1.0-wall_emissivity)>MICRO ) {
      
      
      // Iterate, summing the first moment over the half range for all
      // outgoing directions.  Note that for consistency, sum should
      // add to pi.  This is only possible for boundaries oriented in
      // principle directions.  For arbitrary boundaries, see Modest
      // (2003).  This formulation is used mainly to try and maintain
      // consistency, escpecially for the FVM
      sum = ZERO;
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  sum = sum + max( ZERO, -dot );
	} /* endfor */
      
      // Iterate about each direction and add the absorbed radiation.
      // Add this components energy if this direction is towards the
      // wall (i.e. nx*sx>0 or ny*sy>0)
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  temp = ((1.0-wall_emissivity)/sum);
	  temp *= Uwall.In(v,m,l) * dot;
	  In = In + max( ZERO, temp  );
	  
	  // Note that for the DOM method and the FVM, the term sum=PI
	  // is only consistent for boundaries oriented in the principle
	  // directions.  For arbitrary boundaries, see Modest (2003).

	} /* end for */
      
    } /* endif */
    
    
    //Set the intensity of all outgoing directions, iterate through
    //each direction.
    for (int m=0; m<Uwall.Npolar; m++)
      for (int l=0; l<Uwall.Nazim[m]; l++) {
	
	dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	
	// We are only concerned with directions that are leaving the wall and 
	// entering the domain.  
	if ( dot < ZERO ) Uwall.In(v,m,l) = In;  
	
      }
    
    
  } /* endfor Nband*/

  return (Uwall);

}

/*********************************************************
 * Routine: Gray_Wall_Space_March                        *
 *                                                       *
 * This function sets the wall boundary state.           *
 * Used for Space Marching
 *                                                       *
 *********************************************************/
void Gray_Wall_Space_March(Rte2D_State &Uwall, 
	       const Vector2D &norm_dir, 
	       const double &wall_temperature,
	       const double &wall_emissivity ) {

  double In, dot, temp, sum;
  
 
  // The intensity leaving the wall is the same for all outgoing
  // intensities

  // loop over each band
  for (int v=0; v<Uwall.Nband; v++) {
    
    // for a black wall
    In = wall_emissivity * Ib(wall_temperature);

    // For grey wall.
    if ( fabs(1.0-wall_emissivity)>MICRO ) {
      
      
      // Iterate, summing the first moment over the half range for all
      // outgoing directions.  Note that for consistency, sum should
      // add to pi.  This is only possible for boundaries oriented in
      // principle directions.  For arbitrary boundaries, see Modest
      // (2003).  This formulation is used mainly to try and maintain
      // consistency, escpecially for the FVM
      sum = ZERO;
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  sum = sum + max( ZERO, -dot );
	} /* endfor */
      
      // Iterate about each direction and add the absorbed radiation.
      // Add this components energy if this direction is towards the
      // wall (i.e. nx*sx>0 or ny*sy>0)
      for (int m=0; m<Uwall.Npolar; m++)
	for (int l=0; l<Uwall.Nazim[m]; l++) {
	  dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	  temp = ((1.0-wall_emissivity)/sum);
	  temp *= Uwall.In(v,m,l) * dot;
	  In = In + max( ZERO, temp  );
	  
	  // Note that for the DOM method and the FVM, the term sum=PI
	  // is only consistent for boundaries oriented in the principle
	  // directions.  For arbitrary boundaries, see Modest (2003).

	} /* end for */
      
    } /* endif */
    
    
    //Set the intensity of all outgoing directions, iterate through
    //each direction.
    for (int m=0; m<Uwall.Npolar; m++)
      for (int l=0; l<Uwall.Nazim[m]; l++) {
	
	dot = norm_dir.x * Uwall.mu[m][l] + norm_dir.y * Uwall.eta[m][l];
	
	// We are only concerned with directions that are leaving the wall and 
	// entering the domain.  
	if ( dot < ZERO ) Uwall.In(v,m,l) = In;  
	
      }
    
    
  } /* endfor Nband*/


}




/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 * For reflected incoming rays that do not coincide     *
 * with a discrete direction, the interpolation         *
 * technique of                                         *
 *   J. Liu, H.M. Shang, and Y.S. Cheng, J. Quant.      *
 *   Spectrosc. Radiat. Transfer  66 (2000) 17-33.      *
 * is used.                                             *
 *                                                      *
 ********************************************************/
Rte2D_State Reflect(const Rte2D_State &U, const Vector2D &norm_dir) {

    Vector3D out_dir, dir;
    Vector2D int_dir, in_dir;
    double cos_angle, sin_angle;
    double dcn, dct, dot_prod;
    Rte2D_State Temp(U);  Temp.Zero();
    bool exact_match;
    double arc_len, cos_phi, num, denom;
    int mmm, lll;
 
    // store the closest direction
    // for interpolation when rays do not reflect 
    // exactly onto other rays
    int mm, ll;
    double dotp;
	    

    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */
    
      for (int m=0; m<U.Npolar; m++) 
	for (int l=0; l<U.Nazim[m]; l++) {


	  // set the incoming direction
	  in_dir.x = U.mu[m][l];
	  in_dir.y = U.eta[m][l];
	
	  // don't reflect any directions leaving the reflection plane
	  // and entering the domain
	  // if ( dot(in_dir,norm_dir) < ZERO ) continue;
	  

	  /* Determine the direction cosine's for the frame
	     rotation. */
	  dcn =   in_dir.x*norm_dir.x + in_dir.y*norm_dir.y;
	  dct = - in_dir.x*norm_dir.y + in_dir.y*norm_dir.x;
	  
	  /* Reflect the normal component in the rotated frame. */
	  dcn = - dcn;
	  
	  /* Rotate back to the original Cartesian reference frame. */
	  out_dir.x = dcn*norm_dir.x - dct*norm_dir.y;
	  out_dir.y = dcn*norm_dir.y + dct*norm_dir.x;
	  out_dir.z = U.xi[m][l];


	  // now search for the closest direction
	  exact_match = false;
	  dotp = 0; mm = 0; ll = 0;
	  for (int p=0; p<U.Npolar && !exact_match; p++) 
	    for (int q=0; q<U.Nazim[p] && !exact_match; q++) {

	      // compute the dot product between the direction to check
	      // and the normal direction
	      dot_prod = norm_dir.x * U.mu[p][q] + norm_dir.y * U.eta[p][q];
	      
	      // don't look at any directions leaving the domain and
	      // incident on the reflection plane
	      // if ( dot_prod > ZERO ) continue;
	      
	      
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
		
	      } /* endif */
	      

	    } /* endfor */
	  
	  if (!exact_match)
	    cout << endl
		 << "Warning - Rte2DState.cc::Reflect - Control angle overlap "
		 << "detected.  Add pixelation at the boundaries."
		 << endl
		 << "m = " << setw(4) << m << ", l = " << setw(4) << l << endl
		 << "m' = " << setw(4) << mm << ", l' = " << setw(4) << ll << endl
		 << "mu = " << setw(18) << U.mu[m][l] << ", eta = " << setw(18) << U.eta[m][l] << endl
		 << "mu' = " << setw(18) << U.mu[mm][ll] << ", eta' = " << setw(18) << U.eta[mm][ll];
	  
	  // set the reflected intensity
	  // Remember: directions invarient with wavenumber
	  for (int v=0; v<U.Nband; v++) 
	    Temp.In(v,mm,ll) = U.In(v,m,l);

      
	} /* endfor */ 

  
    /* Return the reflected state. */
    return (Temp);
}


/********************************************************
 * Routine: Reflect_Space_March                         *
 *                                                      *
 * This function sets the reflected solution state      *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 * Used for Space Marching
 *                                                      *
 ********************************************************/
void Reflect_Space_March(Rte2D_State &U, const Vector2D &norm_dir) {

    Vector3D out_dir, dir;
    Vector2D int_dir, in_dir;
    double cos_angle, sin_angle;
    double dcn, dct, dot_prod;
    bool exact_match;
    double arc_len, cos_phi, num, denom;
    int mmm, lll;
 
    // store the closest direction
    // for interpolation when rays do not reflect 
    // exactly onto other rays
    int mm, ll;
    double dotp;
	    

    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */

      for (int m=0; m<U.Npolar; m++) 
	for (int l=0; l<U.Nazim[m]; l++) {

	  // set the incoming direction
	  in_dir.x = U.mu[m][l];
	  in_dir.y = U.eta[m][l];
	
	  // don't reflect any directions leaving the reflection plane
	  // and entering the domain
	  if ( dot(in_dir,norm_dir) < ZERO ) continue;
	  

	  /* Determine the direction cosine's for the frame
	     rotation. */
	  dcn =   in_dir.x*norm_dir.x + in_dir.y*norm_dir.y;
	  dct = - in_dir.x*norm_dir.y + in_dir.y*norm_dir.x;
	  
	  /* Reflect the normal component in the rotated frame. */
	  dcn = - dcn;
	  
	  /* Rotate back to the original Cartesian reference frame. */
	  out_dir.x = dcn*norm_dir.x - dct*norm_dir.y;
	  out_dir.y = dcn*norm_dir.y + dct*norm_dir.x;
	  out_dir.z = U.xi[m][l];


	  // now search for the closest direction
	  exact_match = false;
	  dotp = 0;
	  for (int p=0; p<U.Npolar && !exact_match; p++) 
	    for (int q=0; q<U.Nazim[p] && !exact_match; q++) {

	      // compute the dot product between the direction to check
	      // and the normal direction
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
		
	      } /* endif */
	      

	    } /* endfor */
	  
	  if (!exact_match)
	    cout << endl
		 << "Warning - Rte2DState.cc::Reflect_Space_March - Control angle overlap "
		 << "detected.  Add pixelation at the boundaries.";
	  
	  // set the reflected intensity
	  // Remember: directions invarient with wavenumber
	  for (int v=0; v<U.Nband; v++) 
	    U.In(v,mm,ll) = U.In(v,m,l);

      
	} /* endfor */ 

}


 /**************************************************************************
  ************************* SCATTER FUNCTIONS  *****************************
  **************************************************************************/

/********************************************************
 * Routine: Legendre                                    *
 *                                                      *
 * Compute the nth degree legendre polynomial at x      *
 *                                                      *
 ********************************************************/
double Legendre( const double x, const int n) { 

  // initialize
  double P0, P1, Pn;

  // 0th degree
  P0 = 1;
  if (n == 0) return P0;

  // 1st degree
  P1 = x;
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
double* PhaseFunc( const int type, int &n) {

  // declares
  double* An;

  // setup the expansion coefficients
  switch (type) {
	    
    // for Linear isotropic scattering
    case (RTE2D_SCATTER_ISO):
    default:

      n = 1;
      An = new double[n];
      An[0]  = 1.00000;

     break;

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


    // Forward scattering with the F2 phase function of Kim and Lee (1988)
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


    // Forward scattering with the F3 phase function of Kim and Lee (1988)
    case (RTE2D_SCATTER_F3):
      
      // the degree
      n = 3;

      // create the array and set the constants
      An = new double[n];
      An[0] = 1.00000;
      An[1] = 1.20000;
      An[2] = 0.50000;

      break;


    // Backward scattering with the B1 phase function of Kim and Lee (1988)
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


    // Backward scattering with the B2 phase function of Kim and Lee (1988)
    case (RTE2D_SCATTER_B2):
      
      // the degree
      n = 3;

      // create the array and set the constants
      An = new double[n];
      An[0] =  1.00000;
      An[1] = -1.20000;
      An[2] =  0.50000;

      break;

  }

  // return the array
  return An;

}


 /**************************************************************************
  ************************ EXACT SOLN FUNCTIONS  ***************************
  **************************************************************************/


/********************************************************
 * Routine: exact_cyl                                   *  
 *                                                      *
 * This function is needed for integrating the exact    *
 * solution.  It is passed to                           *
 * SimpsonMultiDim.h:adaptsim().                        *
 ********************************************************/
double exact_cyl(int ndim, double *x, void *params) {

  // declares
  exact_cyl_param P = *(exact_cyl_param *)params;
  double jac;
  double theta;
  double theta_min, theta_max;
  double temp = 0;

  // set out integration parameters
  double u = x[0];
  double psi = x[1];

  // compute integration limits
  double theta_B = atan( ( P.r*cos(psi) + sqrt(ONE-pow(P.r,TWO)*pow(sin(psi),TWO))) 
			 / (P.z+P.c) );
  double theta_C =  atan( ( P.r*cos(psi) + sqrt(ONE-pow(P.r,TWO)*pow(sin(psi),TWO))) 
			  / (P.z-P.c) ) + PI;

  // limit them
  theta_B = max(theta_B, ZERO);
  theta_C = min(theta_C, PI);
  
  // Depending on the term, compute the main integrand and the integration limits
  // First, transform the variable x to the variable limit
  // Transforms the domain 0<=x<=1 to theta_min<=theta<=theta_max
  // for integrating  over variable areas
  // Second, compute the main integratand

  switch (P.term_flag) {
  case 0:
    theta_min = 0;
    theta_max = theta_B;
    theta = (1-u)*theta_min + u*theta_max;
    jac = theta_max - theta_min;
    temp = ( ONE - exp( -P.kappa*(P.z+P.c) / cos(theta) ) ) * sin(theta) * jac;
    break;
  case 1:
    theta_min = theta_B;
    theta_max = theta_C;
    theta = (1-u)*theta_min + u*theta_max;
    jac = theta_max - theta_min;
    // checking for singularities
    if (theta==0 || theta==PI) temp = ONE; // singularity, but limit=1
    else temp = ( ONE - exp( -P.kappa*( P.r*cos(psi) + sqrt(ONE-pow(P.r*sin(psi),TWO) ) ) 
			     / sin(theta) ) );
    temp *= sin(theta) * jac;
    break;
  case 2:
    theta_min = theta_C;
    theta_max = PI;
    theta = (ONE-u)*theta_min + u*theta_max;
    jac = theta_max - theta_min;
    temp = ( ONE - exp( -P.kappa*(P.z-P.c) / cos(theta) ) ) * sin(theta) * jac;
    break;
  }


  // Depending on the coordinate system, 
  // multiply by the apprpriate direction cosine
  switch (P.coord_flag) {

    // computing direction integrated intensity
  case 0: 
    temp *= ONE;
    break;
    // computing r-direction heat flux
  case 1:
    temp *= sin(theta)*cos(psi);
    break;
    //computing z-direction heat flux
  case 2:
    temp *= cos(theta);
    break;
  }

  // return the value
  return temp;

}



/********************************************************
 * Routine: CylindricalEnclosure                        *
 *                                                      *
 * This function computes the exact solution for a      *
 * emitting-absorbing  isothermal medium with isothermal*
 * bounding cold walls.  See Dua and Cheng (1975)       *
 *                                                      *
 ********************************************************/
void CylindricalEnclosure( const double gas_temp,
			   const double c,    // cylinder half-length / Ro
			   const double tau,  // kappa*Ro
			   const double rpos, // r/Ro
			   const double zpos, // z/Ro
			   double &G, 
			   double &qr, 
			   double &qz ){


  simp_function F;        // function struct for integration
  exact_cyl_param params; // function parameters struct for integration
  simp_state S;           // state struct for integration
  simp_params P;          // params struct for integration
  int fevals;             // number of function evaluations
  double val;             // value of integration
  int err;                // error flag returned by adaptsim

  // zero
  G = ZERO;
  qr = ZERO;
  qz = ZERO;
  

  // get the function parameters
  params.r = rpos;
  params.z = zpos;
  params.c = c;
  params.kappa = tau;
  params.term_flag = 0;
  params.coord_flag = 0;

  // setup integration ( allocate memory and set parameters )
  malloc_simp_struc( 2, F, S );
  init_simp_struc( F, S );
  
  // setup function
  F.f = exact_cyl;
  F.params = &params;

  // setup integration parameters
  P.maxevals = 100000000;
  P.tol = 1e-6;


  //----------------------- Total irradiation -----------------------//  
  params.coord_flag = 0;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<3 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = ONE;   // theta - variable
    F.xmin[1] = ZERO; F.xmax[1] = TWO*PI;// psi
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr << "RteState.cc::CylindricalEnclosure - "
	   << "Error integrating G \n";
      cerr << "Error flag: " << err << endl;
      exit (-1);
    }

    // add the contribution
    G += val;

  } /* endfor term_flag */ 

  // normalize by 4 PI * blackbody intensity
  G /= FOUR * PI;


  //---------------------------- R-dir Flux -------------------------//
  params.coord_flag = 1;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<3 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = ONE;   // theta - variable
    F.xmin[1] = ZERO; F.xmax[1] = TWO*PI;// psi
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr << "RteState.cc::CylindricalEnclosure - "
	   << "Error integrating qr \n";
      cerr << "Error flag: " << err << endl;
      exit (-1);
    }

    // add the contribution
    qr += val;

  } /* endfor term_flag */ 

  // normalize by PI * blackbody intensity
  qr /= PI;


  //---------------------------- Z-dir Flux -------------------------//
  params.coord_flag = 2;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<3 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = ONE;   // theta - variable
    F.xmin[1] = ZERO; F.xmax[1] = TWO*PI;// psi
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr << "RteState.cc::CylindricalEnclosure - "
	   << "Error integrating qr \n";
      cerr << "Error flag: " << err << endl;
      exit (-1);
    }

    // add the contribution
    qz += val;

  } /* endfor term_flag */ 

  // normalize by PI * blackbody intensity
  qz /= PI;

}



/********************************************************
 * Routine: exact_rect                                  *  
 *                                                      *
 * This function is needed for integrating the exact    *
 * solution.  It is passed to                           *
 * SimpsonMultiDim.h:adaptsim().                        *
 ********************************************************/
double exact_rect(int ndim, double *x, void *params) {

  // declares
  exact_rect_param P = *(exact_rect_param *)params;
  double jac;
  double psi, psi_min, psi_max;
  double s1, s1_star, s1_til;
  double d, x_star, y_star;

  // set out integration parameters
  double theta = x[0];
  double u = x[1];

  // compute integration limits
  double psi_A = atan( (P.y-P.b2)/(P.x-P.a1) );
  double psi_B = atan( (P.y-P.b1)/(P.x-P.a1) );
  double psi_C = atan( (P.y-P.b1)/(P.x-P.a2) );
  double psi_D = atan( (P.y-P.b2)/(P.x-P.a2) );


  // Depending on the term, compute the main integrand 
  // and the integration limits
  // 1 -> Transform the variable u to the variable limit psi  
  //      ( <=u<=1 to theta_min<=theta<=theta_max)
  // 2 -> rotate the coordinate frame
  // 3 -> compute the main integrand
  double temp = 0;
  switch (P.term_flag) {

  case 0:
    // variable transformation
    psi_max = psi_B;
    psi_min = psi_A;
    psi = (1-u)*psi_min + u*psi_max;
    jac = psi_max - psi_min;

    // field point coordinate in rotated frame
    s1 =   P.x*cos(psi) + P.y*sin(psi);

    // distance to west wall along ray
    d = ( P.x - P.a1 ) / ( sin(theta) * cos(psi) );

    // ray emission wall point coordinate in rotated frame
    x_star = P.a1;
    y_star = P.y - d*sin(theta)*sin(psi); 
    s1_star = x_star*cos(psi) + y_star*sin(psi);
    break;

  case 1:
    // variable transformation
    psi_max = psi_D+PI;
    psi_min = psi_C+PI;
    psi = (1-u)*psi_min + u*psi_max;
    jac = psi_max - psi_min;

    // rotate our coordinate frame
    s1 = P.x*cos(psi) + P.y*sin(psi);

    // distance to east wall along ray
    d = ( P.x - P.a2 ) / ( sin(theta) * cos(psi) );

    // ray emission wall point coordinate in rotated frame
    x_star = P.a2;
    y_star = P.y - d*sin(theta)*sin(psi); 
    s1_star = x_star*cos(psi) + y_star*sin(psi);
    break;

  case 2:
    // variable transformation
    psi_max = psi_C+PI;
    psi_min = psi_B;
    psi = (1-u)*psi_min + u*psi_max;
    jac = psi_max - psi_min;

    // rotate our coordinate frame
    s1 = P.x*cos(psi) + P.y*sin(psi);

    // distance to south wall along ray
    d = ( P.y - P.b1 ) / ( sin(theta) * sin(psi) );

    // ray emission wall point coordinate in rotated frame
    y_star = P.b1;
    x_star = P.x - d*sin(theta)*cos(psi); 
    s1_star = x_star*cos(psi) + y_star*sin(psi);
    break;

  case 3:
    // variable transformation
    psi_max = psi_A+2*PI;
    psi_min = psi_D+PI;
    psi = (1-u)*psi_min + u*psi_max;
    jac = psi_max - psi_min;

    // rotate our coordinate frame
    s1 = P.x*cos(psi) + P.y*sin(psi);

    // distance to north wall along ray
    d = ( P.y - P.b2 ) / ( sin(theta) * sin(psi) );

    // ray emission wall point coordinate in rotated frame
    y_star = P.b2;
    x_star = P.x - d*sin(theta)*cos(psi); 
    s1_star = x_star*cos(psi) + y_star*sin(psi);
      
    break;
  }

  // some integration constants
  s1_til = (s1 - s1_star) / sin(theta);

  // compute integrand
  // first, check if the ray actually originated from a point on the wall
  if ( x_star>=P.a1 && x_star<=P.a2 && y_star>=P.b1 && y_star<=P.b2) {
    temp =  exp( P.kappa * (d-s1_til) );
    temp -= exp( -P.kappa*s1_til );
  } else {
    temp = 0;
  }
  
  // multiply by integration parameters
  temp *= sin(theta) * jac;


  // Depending on the coordinate system, 
  // multiply by the apprpriate direction cosine
  switch (P.coord_flag) {

    // computing direction integrated intensity
  case 0: 
    temp *= ONE;
    break;
    // computing x-direction heat flux
  case 1:
    temp *= sin(theta)*cos(psi);
    break;
    //computing y-direction heat flux
  case 2:
    temp *= sin(theta)*sin(psi);
    break;
  }

  // return the value
  return temp;

}



/********************************************************
 * Routine: RectangularEnclosure                        *
 *                                                      *
 * This function computes the exact solution for a      *
 * emitting-absorbing  isothermal medium with isothermal*
 * bounding cold walls.  See Cheng (1972)       *
 *                                                      *
 ********************************************************/
void RectangularEnclosure( const double gas_temp,
			   const double kappa,  // absorbsion coeff
			   const double left,   // west wall location
			   const double right,  // east wall location
			   const double bot,    // south wall location
			   const double top,    // north wall location
			   const double xpos,   // x position
			   const double ypos,   // y position
			   double &G, 
			   double &qx, 
			   double &qy ){


  simp_function F;        // function struct for integration
  exact_rect_param params; // function parameters struct for integration
  simp_state S;           // state struct for integration
  simp_params P;          // params struct for integration
  int fevals;             // number of function evaluations
  double val;             // value of integration
  int err;                // error flag returned by adaptsim

  // zero
  G = ZERO;
  qx = ZERO;
  qy = ZERO;
  

  // get the function parameters
  params.x = xpos;
  params.y = ypos;
  params.a1 = left;
  params.a2 = right;
  params.b1 = bot;
  params.b2 = top;
  params.kappa = kappa;
  params.term_flag = 0;
  params.coord_flag = 0;

  // setup integration ( allocate memory and set parameters )
  malloc_simp_struc( 2, F, S );
  init_simp_struc( F, S );
  
  // setup function
  F.f = exact_rect;
  F.params = &params;

  // setup integration parameters
  P.maxevals = 100000;
  P.tol = 1e-6;


  //----------------------- Total irradiation -----------------------//  
  params.coord_flag = 0;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<4 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = PI;   // theta
    F.xmin[1] = ZERO; F.xmax[1] = ONE;  // psi - variable
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr.precision(4);
      cerr << "RteState.cc::RectangularEnclosure - "
	   << "Error integrating G \n";
      cerr << "Error flag: " << err << endl;
      cerr << "Value:      " << val << endl;
      cerr << "Point:      (" << xpos << ", " << ypos << ")" << endl;
      exit (-1);
    }

    // add the contribution
    G += val;

  } /* endfor term_flag */ 

  // normalize by 4 PI * blackbody intensity
  G /= FOUR * PI;


  //---------------------------- x-dir Flux -------------------------//
  params.coord_flag = 1;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<4 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = PI;   // theta
    F.xmin[1] = ZERO; F.xmax[1] = ONE;  // psi - variable
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr.precision(4);
      cerr << "RteState.cc::RectangularEnclosure - "
	   << "Error integrating qx \n";
      cerr << "Error flag: " << err << endl;
      cerr << "Value:      " << val << endl;
      cerr << "Point:      (" << xpos << ", " << ypos << ")" << endl;
      exit (-1);
    }

    // add the contribution
    qx += val;

  } /* endfor term_flag */ 

  // normalize by PI * blackbody intensity
  qx /= PI;


  //---------------------------- y-dir Flux -------------------------//
  params.coord_flag = 2;
  
  // Use an adaptive simpsons quadrature rule 
  // for multidimensional integration

  // loop over each term
  for(int i=0 ; i<4 ; i++) {

    // initialize before integration
    init_simp_struc( F, S );
    
    // set new integration limits
    F.xmin[0] = ZERO; F.xmax[0] = PI;   // theta
    F.xmin[1] = ZERO; F.xmax[1] = ONE;  // psi - variable
  
    // set the term flag
    params.term_flag = i;
	  
    // compute the integrated phase function
    err = adaptsim( F, S, P, fevals, val );
    if (err) { 
      cerr.precision(4);
      cerr << "RteState.cc::RecatngularEnclosure - "
	   << "Error integrating qy \n";
      cerr << "Error flag: " << err << endl;
      cerr << "Value:      " << val << endl;
      cerr << "Point:      (" << xpos << ", " << ypos << ")" << endl;
      exit (-1);
    }

    // add the contribution
    qy += val;

  } /* endfor term_flag */ 

  // normalize by PI * blackbody intensity
  qy /= PI;

}

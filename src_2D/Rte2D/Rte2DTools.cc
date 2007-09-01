/******************* Rte2DTools.cc ************************************

  Header file for Various functions for Rte2D that do not really 
  belong anywhere.

***********************************************************************/
#include "Rte2DTools.h"

/**************************************************************************
 **************************** SETUP FUNCTIONS  ****************************
 **************************************************************************/

/*******************************************************
 * Main function to setup the Rte2D_State static       *
 * parameters                                          *
 *******************************************************/  
void SetupStateStatic( const Rte2D_Input_Parameters &IP ) 
{

  // deallocate static vars first, just to be safe
  Rte2D_State::DeallocateStatic();

  //------------------------------------------------
  // Absorbsion 
  //------------------------------------------------
  // set the absorption type flag
  Rte2D_State::Absorb_Type = IP.i_AbsorptionModel;

  // GRAY
  if (IP.i_AbsorptionModel == RTE2D_ABSORB_GRAY) {
    Rte2D_State::Nband = 1;

  // SNBCK
  } else if (IP.i_AbsorptionModel == RTE2D_ABSORB_SNBCK) {
    Rte2D_State::AllocateSNBCK();
    Rte2D_State::SNBCKdata->Setup(IP.SNBCK_IP, IP.CFFC_Path);
    Rte2D_State::Nband = Rte2D_State::SNBCKdata->NumVar();

  // ERROR
  } else {
    cerr << "Rte2D_State::SetupState - Invalid flag for gas type\n";
    exit(-1);
  } // endif


  //------------------------------------------------
  // Directions, Phase function 
  //------------------------------------------------
  // FVM
  if (IP.i_RTE_Solver == RTE2D_SOLVER_FVM) {
    Rte2D_State::SetDirsFVM(IP.Number_of_Angles_Mdir, 
			    IP.Number_of_Angles_Ldir, 
			    IP.Axisymmetric );
    Rte2D_State::SetupPhaseFVM( IP.i_ScatteringFunc );

  // DOM
  } else if (IP.i_RTE_Solver == RTE2D_SOLVER_DOM) {
    Rte2D_State::SetDirsDOM(IP.i_DOM_Quadrature,
			    IP.Axisymmetric, 
			    IP.CFFC_Path );
    Rte2D_State::SetupPhaseDOM( IP.i_ScatteringFunc );

  // ERROR
  } else {
    cerr << "SetupStateStatic() - Invalid flag for RTE solver\n";
    exit(1);
  } // endif


}

/*******************************************************
 * Depending on the input parameters, initilize the    *
 * non-solution parameters of state U accordingly.     *
 *******************************************************/  
void SetInitialValues( Rte2D_State &U, 
		       const Rte2D_Input_Parameters &IP ) 
{

  //------------------------------------------------
  // Absorbsion coefficient, Blackbody intensity 
  //------------------------------------------------
  // Use SNBCK
  if (IP.i_AbsorptionModel == RTE2D_ABSORB_SNBCK) {
    U.SNBCKdata->CalculateAbsorb( IP.Pressure/PRESSURE_STDATM, //[atm]
				  IP.Temperature,              //[K]
				  IP.xco,
				  IP.xh2o,
				  IP.xco2,
				  IP.xo2,
				  IP.fsoot,
				  U.kappa );
    U.SNBCKdata->CalculatePlanck( IP.Temperature, U.Ib );

  // Use Gray Gas (ie. constant)
  } else if (IP.i_AbsorptionModel == RTE2D_ABSORB_GRAY) {
    U.SetAbsorption( IP.AbsorptionCoef );
    U.SetBlackbody( BlackBody(IP.Temperature) );

  // error
  } else{
    cerr << "SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  //------------------------------------------------
  // Scattering coefficient 
  //------------------------------------------------
  // scattering coefficient always assumed gray
  U.SetScattering( IP.ScatteringCoef );
  
}


 /**************************************************************************
  ************************ EXACT SOLN FUNCTIONS  ***************************
  **************************************************************************/


/********************************************************
 * Routine: exact_cyl                                   *  
 *                                                      *
 * This function is needed for integrating the exact    *
 * solution.  In CylindricalEnclosure, it is passed to  *
 * SimpsonMultiDim.h:adaptsim() and is integrated.      *
 *                                                      *
 * See Dua and Cheng (1975)                             *
 *                                                      *
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
  
  //------------------------------------------------
  // Depending on the term, compute the main integrand and the integration limits.
  // Each parameter that we are computing, G, qr, qz, consists of three integral terms.
  // They are defined on one of the three intervals:
  //    0<=theta<theta_B,  theta_B<=theta<theta_C,  theta_C<=theta<PI
  // P.term_flag is specifies which term we are evaluating.

  // Perform the following
  // 1. Transform the variable x on the domain 0<=x<=1 to theta_min<=theta<=theta_max.  
  //    This is required for integrating  over variable areas.
  // 2. Compute the main integrand term
  //------------------------------------------------
  switch (P.term_flag) {
    // first integral
  case 0: 
    theta_min = 0;
    theta_max = theta_B;
    theta = (1-u)*theta_min + u*theta_max;
    jac = theta_max - theta_min;
    temp = ( ONE - exp( -P.kappa*(P.z+P.c) / cos(theta) ) ) * sin(theta) * jac;
    break;
    // second integral
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
    // third integral
  case 2:
    theta_min = theta_C;
    theta_max = PI;
    theta = (ONE-u)*theta_min + u*theta_max;
    jac = theta_max - theta_min;
    temp = ( ONE - exp( -P.kappa*(P.z-P.c) / cos(theta) ) ) * sin(theta) * jac;
    break;
  } // endswitch


  //------------------------------------------------
  // Depending on the coordinate system, (denoted by P.coord_flag)
  // multiply by the apprpriate direction cosine
  //------------------------------------------------
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
  } // endswitch

  // return the value
  return temp;

}



/********************************************************
 * Routine: CylindricalEnclosure                        *
 *                                                      *
 * This function computes the exact solution for a      *
 * emitting-absorbing  isothermal medium within an      *
 * black, isothermal enclosure with cold walls.  An     *
 * adaptive simpsons quadrature rule for                *
 * multidimensional integration is used.                *
 *                                                      *
 * See Dua and Cheng (1975)                             *
 *                                                      *
 ********************************************************/
  // Use an adaptive simpsons quadrature rule for multidimensional 
  // integration.
void CylindricalEnclosure( const double gas_temp,
			   const double c,    // cylinder half-length / Ro
			   const double tau,  // kappa*Ro
			   const double rpos, // r/Ro
			   const double zpos, // z/Ro
			   double &G, 
			   double &qr, 
			   double &qz ){

  // declares
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

  // Setup integration: allocate memory and set parameters.
  malloc_simp_struc( 2, F, S );
  init_simp_struc( F, S );
  
  // setup function
  F.f = exact_cyl;
  F.params = &params;

  // setup integration parameters
  P.maxevals = 100000000;
  P.tol = 1e-6;


  //------------------------------------------------
  // Total irradiation   
  //------------------------------------------------
  params.coord_flag = 0;
  
  //
  // loop over each term
  //
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

  } // endfor term_flag

  // normalize by 4 PI * blackbody intensity
  G /= FOUR * PI;


  //------------------------------------------------
  // R-dir Flux 
  //------------------------------------------------
  params.coord_flag = 1;
  
  //
  // loop over each term
  //
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

  } // endfor term_flag 

  // normalize by PI * blackbody intensity
  qr /= PI;


  //------------------------------------------------
  // Z-dir Flux 
  //------------------------------------------------
  params.coord_flag = 2;
  
  //
  // loop over each term
  //
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

  } // endfor term_flag 

  // normalize by PI * blackbody intensity
  qz /= PI;

  // free memory
  free_simp_struc( F, S );

}



/********************************************************
 * Routine: exact_rect                                  *  
 *                                                      *
 * This function is needed for integrating the exact    *
 * solution.  In RectangularEnclosure, it is passed to  *
 * SimpsonMultiDim.h:adaptsim() and integrated.         *
 *                                                      *
 * See Cheng (1972)                                     *
 *                                                      *
 ********************************************************/
double exact_rect(int ndim, double *x, void *params) {

  // declares
  exact_rect_param P = *(exact_rect_param *)params;
  double jac;
  double psi, psi_min, psi_max;
  double s1, s1_star, s1_til;
  double d, x_star, y_star;
  double temp = ZERO;

  // set out integration parameters
  double theta = x[0];
  double u = x[1];

  // compute integration limits
  double psi_A = atan( (P.y-P.b2)/(P.x-P.a1) );
  double psi_B = atan( (P.y-P.b1)/(P.x-P.a1) );
  double psi_C = atan( (P.y-P.b1)/(P.x-P.a2) );
  double psi_D = atan( (P.y-P.b2)/(P.x-P.a2) );


  //------------------------------------------------
  // Depending on the term, compute the main integrand and the integration limits.
  // Each parameter that we are computing, G, qx, qy, consists of four integral terms.
  // They are defined on one of the three intervals:
  //    0<=theta<theta_B,  theta_B<=theta<theta_C,  theta_C<=theta<PI
  // P.term_flag is specifies which term we are evaluating.

  // Perform the following
  // 1. Transform the variable x on the domain 0<=x<=1 to theta_min<=theta<=theta_max.  
  //    This is required for integrating  over variable areas.
  // 2. Compute the main integrand term
  //------------------------------------------------
  switch (P.term_flag) {

  //
  // first term
  //
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

  //
  // second term
  //
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

  //
  // third term
  //
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

  //
  // fourth term
  //
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
  } // endswitch


  //
  // compute integrand
  //
  // first, check if the ray actually originated from a point on the wall
  if ( x_star>=P.a1 && x_star<=P.a2 && y_star>=P.b1 && y_star<=P.b2) {
    s1_til = (s1 - s1_star) / sin(theta);
    temp =  exp( P.kappa * (d-s1_til) );
    temp -= exp( -P.kappa*s1_til );
  } else {
    temp = 0;
  }
  
  // multiply by integration parameters
  temp *= sin(theta) * jac;


  //------------------------------------------------
  // Depending on the coordinate system, 
  // multiply by the apprpriate direction cosine
  //------------------------------------------------
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
 * emitting-absorbing  isothermal medium within an      *
 * black, isothermal enclosure with cold walls.  An     *
 * adaptive simpsons quadrature rule for                *
 * multidimensional integration is used.                *
 *                                                      *
 * See Cheng (1972)                                     *
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

  // declines
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


  //------------------------------------------------
  // Total irradiation 
  //------------------------------------------------
  params.coord_flag = 0;

  //
  // loop over each term
  //
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


  //------------------------------------------------
  // x-dir Flux 
  //------------------------------------------------
  params.coord_flag = 1;
  
  //
  // loop over each term
  //
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


  //------------------------------------------------
  // y-dir Flux 
  //------------------------------------------------
  params.coord_flag = 2;
  
  //
  // loop over each term
  //
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

  // free memory
  free_simp_struc( F, S );

}
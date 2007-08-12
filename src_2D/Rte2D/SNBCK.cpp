/*******************************************************************
  File: SNBCK.cpp

  Description:  This file defines the SNBCK class member functions for 
                the computation of the absorbsion coefficient using
                the SNBCK.  

  Author:  Marc R.J. Charest

  Date:    Dec 25th, 2006
*******************************************************************/
#include "SNBCK.h"



/*********************************************************************
 ********************* EM2C CLASS MEMBER FUNCTIONS *******************
 *********************************************************************/
 
/*********************************************************************
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/


/*********************************************************************
 * EM2C :: LoadParams                                                *
 *                                                                   *
 * Read SNB model parameters, see Soufiani and Taine (1997).         *
 * SNB parameters k (in cm-1 atm-1) and (1/delta) (in cm) are given  *
 * in the files SNBH2O, SNBCO2 and SNBCO for H2O, CO2 and CO         *
 * respectively with a 25 cm-1 spectral resolution. For H2O, the     *
 * parameters are given for bands centred at a wave number ranging   *
 * continuously from 150 cm-1 to 9300 cm-1 (367 bands). For the      *
 * other species, the parameters are not given in the transparency   *
 * regions. For CO2, the centers of the contributing bands range     *
 * from 450 to 1200 cm-1 (31 bands), from 1950 to 2450 cm-1          *
 * (21 bands), from 3300 to 3800 cm-1 (21 bands) and from 4700 to    *
 * 5250 cm-1 (23 bands). The total number of CO2 bands is 96. For CO,*
 * the centers of the contributing bands range from 1750 to 2325 cm-1*
 * (24 bands) and from 3775 to 4350 cm-1 (24 bands). The total number*
 * of CO bands is 48.  The parameters are given for the temperature  *
 * range (300, 2900K) with a constant temperature step T=200K       *
 * (14 temperatures).  Lastly, the file SNBWN contains the spectral  *
 * position of the band center and the band width for the 367 bands. *
 *                                                                   *
 * Inputs: PATH - the path to the folder containing the datafiles.   *
 *********************************************************************/
void EM2C :: LoadParams(const char *PATH) {

  // declares
  ifstream in;
  string str;


  //------------------------------------------------
  // CO PARAMS 
  //------------------------------------------------

  strcpy(FileSNBCO,PATH);
  strcat(FileSNBCO,"/SNBCO");

  // open file 
  in.open(FileSNBCO, ios::in);

  // read parameters and build a natural spline
  // Note: Need to convert FORTRAN "D" notation to CPP "e" notation
  // for exponentials

  // SNB parameter k
  for (int i=0; i<NCO; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      kCO[i][j] = FortranToCppD(str);
    }
  }

  // SNB parameter (1/delta)
  for (int i=0; i<NCO; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      dCO[i][j] = FortranToCppD(str);
    }
  }

  // close file
  in.close();


  //------------------------------------------------
  // H2O PARAMS
  //------------------------------------------------

  strcpy(FileSNBH2O,PATH);
  strcat(FileSNBH2O,"/SNBH2O");

  // open file 
  in.open(FileSNBH2O, ios::in);

  // read parameters and build a natural spline
  // Note: Need to convert FORTRAN "D" notation to CPP "e" notation
  // for exponentials

  // SNB parameter k
  for (int i=0; i<NH2O; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      kH2O[i][j] = FortranToCppD(str);
    }
  }

  // SNB parameter (1/delta)
  for (int i=0; i<NH2O; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      dH2O[i][j] = FortranToCppD(str);
    }
  }

  // close file
  in.close();




  //------------------------------------------------
  // CO2 PARAMS 
  //------------------------------------------------

  strcpy(FileSNBCO2,PATH);
  strcat(FileSNBCO2,"/SNBCO2");

  // open file 
  in.open(FileSNBCO2, ios::in);

  // read parameters and build a natural spline
  // Note: Need to convert FORTRAN "D" notation to CPP "e" notation
  // for exponentials

  // SNB parameter k
  for (int i=0; i<NCO2; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      kCO2[i][j] = FortranToCppD(str);
    }
  }

  // SNB parameter (1/delta)
  for (int i=0; i<NCO2; i++ ) {
    for (int j=0; j<Npoints; j++ ){
      in >> setw(13) >> str;
      dCO2[i][j] = FortranToCppD(str);
    }
  }

  // close file
  in.close();


  //------------------------------------------------
  // WAVENUMBER CENTERS 
  //------------------------------------------------

  strcpy(FileSNBWN,PATH);
  strcat(FileSNBWN,"/SNBWN");

  // open file 
  in.open(FileSNBWN, ios::in);

  // read parameters
  // wavenumber center, bandwidth
  for (int i=0; i<Nbands; i++ ) {
    in >> WaveNo[i] >> BandWidth[i];
  }

  // close file
  in.close();


  // Determin spectral indexes
  FindIndex();


}


/*********************************************************************
 * EM2C :: FindIndex                                                 *
 *                                                                   *
 * Determine Spectral Indexes relating the overall band index to the *
 * index for the data arrays kXX[][] and dXX[][] where XX is the     *
 * species name. This is taken from the original EM2C fortran        *
 * program provided by Soufiani and Taine (1997).  For consistency 1 *
 * is subtracted  from the indexes because Fortran arrays index      *
 * starting from 1.                                                  *
 *********************************************************************/
void EM2C :: FindIndex() {

  // loop over the number of bands
  for (int i=0; i<Nbands; i++) {

    // CO2
    if (WaveNo[i]>=450.0 && WaveNo[i]<=1200)
      iCO2[i]=int((WaveNo[i]-450.0)/25.0)+1-1;
    else if (WaveNo[i]>=1950.0 && WaveNo[i]<=2450.0)
      iCO2[i]=int((WaveNo[i]-1950.0)/25.0)+32-1;
    else if (WaveNo[i]>=3300.0 && WaveNo[i]<=3800.0)
      iCO2[i]=int((WaveNo[i]-3300.0)/25.0)+53-1;
    else if (WaveNo[i]>=4700.0 && WaveNo[i]<=5250.0)
      iCO2[i]=int((WaveNo[i]-4700.0)/25.0)+74-1;
    else 
      iCO2[i]=-9;


    // H2O
    if (WaveNo[i]>=150.0 && WaveNo[i]<=9300.0)
      iH2O[i]=int((WaveNo[i]-150.0)/25.0)+1-1;
    else
      iH2O[i]=-9;

    // CO
    if (WaveNo[i]>=1750.0 && WaveNo[i]<=2325.0)
      iCO[i]=int((WaveNo[i]-1750.0)/25.0)+1-1;
    else if (WaveNo[i]>=3775.0 && WaveNo[i]<=4350.0)
      iCO[i]=int((WaveNo[i]-3775.0)/25.0)+25-1;
    else 
      iCO[i]=-9;

    // set the flag, 1 for active band, 0 for inactive
    liCO[i]  = ( iCO[i]>=0 ? 1 : 0 );
    liCO2[i] = ( iCO2[i]>=0 ? 1 : 0 );
    liH2O[i] = ( iH2O[i]>=0 ? 1 : 0 );
    liMix[i] = ( liCO[i]+liCO2[i]+liH2O[i]>0 ? 1 : 0 );

  } /* end for */




}

/*********************************************************************
 * EM2C :: FindIndex                                                 *
 *                                                                   *
 * Calculate piecewise linear interpolation coefficients.            *
 * This is taken from the original EM2C fortran program provided     *
 * by Soufiani and Taine (1997).  Note that parameters are given for *
 * the temperature range (300,2900K) with a constant temperature     *
 * step of 200 K (i.e. 14 temperatures). For consistency 1           *
 * is subtracted  from the indexes because Fortran arrays index      *
 * starting from 1.                                                  *
 *                                                                   *
 * Inputs:  T - the temperature [K]                                  *
 * Outputs: slope - the slope of the linear interpolation line.      *
 *          index - the index of the left bracketing data point      *
 *********************************************************************/
void EM2C :: Interp( const double T, double &slope, int &index) {

  // if inside bracket, then linearly interpolate
  if(T>300.0 && T<2900.0) {
    slope=(T-300.0)/200.0;
    index=int(slope+MICRO);
    slope=slope-index;
    //index=index+1;

  // to the right of the bracket
  } else if (T>=2900.0) {
    slope=ONE;
    index=12; // last index is 13 -> then one before it is 12

  // to the left of the bracket
  } else {
    slope=ZERO;
    index=0; // fisrt index is 0
  }/* endif */

}
  
/*********************************************************************
 * EM2C :: ComputeRefSNB                                             *
 *                                                                   *
 * Compute SNB model parameters (S and B) for a specified mixture    *
 * (mole fractions of CO, CO2, H2O), temperature, and pressure.      *
 * These parameters are normalized for a unit mole fraction          *
 * and pressure assuming B is independant of mixture composition and *
 * pressure.  See Liu and Smallwood (2004).                          *
 * Piecewise linear interpolants are used.                           *
 *********************************************************************/
void EM2C :: ComputeRefSNB( const double p,       // pressure [atm]
			    const double T,       // temperature [K]
			    const double xCO,     // mole fraction oh CO
			    const double xH2O,    // mole fraction oh H2O
			    const double xCO2,    // mole fraction oh CO2
			    const double xO2 )    // mole fraction oh O2
{
   
  // declares
  int index;
  double slope;
  double gamma, k, delta_inv, beta;

  // check if temperature is within range
  if (T<Tmin || T>Tmax) {
    cerr << "EM2C::ComputeSNB(): Temperature out of valid range, "
	 << Tmin << " < T < " << Tmax << " where T = "
	 << T << " [K].\n";
  }
  

  // compute interpolation coefficients using linear piecewise interpolation
  Interp( T, slope, index);

  // some coeffients
  double T296=296.0/T;
  double T273=273.0/T;
  double T900=900.0/T;
  double xN2 = ONE-xCO-xCO2-xH2O;


  // loop over the number of bands
  for (int i=0; i<Nbands; i++) {

    // zero 
    S_CO[i]  = ZERO; B_CO[i]  = ZERO;
    S_CO2[i] = ZERO; B_CO2[i] = ZERO;
    S_H2O[i] = ZERO; B_H2O[i] = ZERO;

    
    //-------------------- CO -----------------------//
    if (liCO[i]) {

      // Eq. (2) in Soufiani and Taine (1997)
      gamma = 0.075*xCO2*pow(T296,0.60) +
	      0.120*xH2O*pow(T296,0.82) +
              0.060*pow(T296,0.7)*(ONE-xCO2-xH2O);
      gamma *= p;

      // k
      k = kCO[iCO[i]][index] + 
          slope*(kCO[iCO[i]][index+1]-kCO[iCO[i]][index]);

      // d
      delta_inv = dCO[iCO[i]][index] + 
                  slope*(dCO[iCO[i]][index+1]-dCO[iCO[i]][index]);

      // beta
      beta = TWO*gamma*delta_inv;
      
      // S, B
      S_CO[i] = k;
      B_CO[i] = TWO*beta/PI;
      
//       // check for infinite pressure broadening
//        if (B_CO[i]>GRAY_GAS_CUTOFF) 
// 	cerr << "Gray gas (co) detected at " << WaveNo[i] << " cm^-1" << endl;

    }/* endif CO */

    //-------------------- CO2 -----------------------//
    if (liCO2[i]) {

      // Eq. (2) in Soufiani and Taine (1997)
      gamma = 0.070*xCO2 + 0.058*(ONE-xCO2-xH2O) + 0.1*xH2O;
      gamma *= p*pow(T296,0.7);

      // k
      k = kCO2[iCO2[i]][index] + 
          slope*(kCO2[iCO2[i]][index+1]-kCO2[iCO2[i]][index]);

      // d
      delta_inv = dCO2[iCO2[i]][index] + 
                  slope*(dCO2[iCO2[i]][index+1]-dCO2[iCO2[i]][index]);

      // beta
      beta = TWO*gamma*delta_inv;
      
      // S, B
      S_CO2[i] = k;
      B_CO2[i] = TWO*beta/PI;

//       // check for infinite pressure broadening
//        if (B_CO2[i]>GRAY_GAS_CUTOFF) 
// 	cerr << "Gray gas (co2) detected at " << WaveNo[i] << " cm^-1" << endl;
     
    }/* endif CO2 */


    //-------------------- H2O -----------------------//
    if (liH2O[i]) {

      // Eq. (2) in Soufiani and Taine (1997)
      gamma = 0.462*xH2O*T296 +
	      sqrt(T296)*( 0.079*(ONE-xCO2-xO2) + 
			   0.106*xCO2 + 0.036*xO2 );
      gamma *= p;

      // k
      k = kH2O[iH2O[i]][index] + 
	  slope*(kH2O[iH2O[i]][index+1]-kH2O[iH2O[i]][index]);

      // d
      delta_inv = dH2O[iH2O[i]][index] + 
	          slope*(dH2O[iH2O[i]][index+1]-dH2O[iH2O[i]][index]);

      // beta
      beta = TWO*gamma*delta_inv;
      
      // S, B
      S_H2O[i] = k;
      B_H2O[i] = TWO*beta/PI;
      
//       // check for infinite pressure broadening
//       if (B_H2O[i]>GRAY_GAS_CUTOFF) 
// 	cerr << "Gray gas (h2o) detected at " << WaveNo[i] << " cm^-1" << endl;

    }/* endif H2O */

    //-----------------------------------------------//

  } /* end for */

}

/*********************************************************************
 * EM2C :: ComputeRefSNB                                             *
 *                                                                   *
 * Compute SNB model parameters (S and B) for a specified mixture    *
 * (mole fractions of CO, CO2, H2O), temperature, and pressure.      *
 *********************************************************************/
void EM2C :: ComputeSNB( const double p,       // pressure [atm]
			 const double T,       // temperature [K]
			 const double xCO,     // mole fraction oh CO
			 const double xH2O,    // mole fraction oh H2O
			 const double xCO2,    // mole fraction oh CO2
			 const double xO2 ) {  // mole fraction oh O2
  bool flag;

  // compute normalized parameters
  ComputeRefSNB( p, T, xCO, xH2O, xCO2, xO2 );

  // now multiply by species concentration and pressure
  for (int i=0; i<Nbands; i++) {

    // ------------- Individual Species ------------//
    S_CO[i]  *=  xCO*p;
    S_CO2[i] *= xCO2*p;
    S_H2O[i] *= xH2O*p;


    // ---------- Optically Thin Mixture -----------//
    // Treat mixture overlapping bands are treated using the 
    // optically thin approximation proposed by Liu et al. (2001) 
    // by treating the overlapping bands to be an approximate Malkmus band.
    // ---------- Optically Thick Mixture -----------//
    // Treat mixture overlapping bands are treated using the 
    // optically thick approximation proposed by Liu et al. (2001) 
    // by treating the overlapping bands to be an approximate Malkmus band.
     S_Thin[i] = ZERO;  B_Thin[i] = ZERO;
    S_Thick[i] = ZERO; B_Thick[i] = ZERO;
    flag = false;
    if (liH2O[i] && xH2O>MICRO) {
      S_Thin[i] += S_H2O[i];
      B_Thin[i] += pow(S_H2O[i],2)/B_H2O[i];
      S_Thick[i] += S_H2O[i];
      B_Thick[i] += sqrt(S_H2O[i]*B_H2O[i]);
      flag = true;
    } 
    if (liCO2[i] && xCO2>MICRO) {
      S_Thin[i] += S_CO2[i];
      B_Thin[i] += pow(S_CO2[i],2)/B_CO2[i];
      S_Thick[i] += S_CO2[i];
      B_Thick[i] += sqrt(S_CO2[i]*B_CO2[i]);
      flag = true;
    }
    if (liCO[i] && xCO>MICRO) {
      S_Thin[i] += S_CO[i];
      B_Thin[i] += pow(S_CO[i],2)/B_CO[i];
      S_Thick[i] += S_CO[i];
      B_Thick[i] += sqrt(S_CO[i]*B_CO[i]);
      flag = true;
    }
   
    // compute SNB parameter, check for div by zero
    if ( flag ) {
      B_Thin[i] = pow(S_Thin[i],2)/B_Thin[i]; 
      B_Thick[i] = pow(B_Thick[i],2)/S_Thick[i];
    }


  } /* end for */

}


/*********************************************************************
 * EM2C :: Transmissivity                                            *
 *                                                                   *
 * Computes the transmissivity over an isothermal and homogeneous    *
 * path assuming a Malkmus gas. Eq (4) in:                           *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 *                                                                   *
 *      L - path length in cm                                        *
 *********************************************************************/
double EM2C :: Transmissivity( const double L, const int i ){

  // transparent case
  double tau = ONE;
  if (liCO2[i]) tau *= TransmissivitySNB( B_CO2[i], S_CO2[i], L );
  if (liH2O[i]) tau *= TransmissivitySNB( B_H2O[i], S_H2O[i], L );
  if (liCO[i])  tau *= TransmissivitySNB(  B_CO[i],  S_CO[i], L );
  return tau;

}

/*********************************************************************
 ********************* SNBCK CLASS MEMBER FUNCTIONS ******************
 *********************************************************************/

/*********************************************************************
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/


/*********************************************************************
 * Constructor that initializes without values                       *
 *********************************************************************/
SNBCK :: SNBCK() :  SNB(), Nquad(0), Nbands(0),  
		    g(NULL), w(NULL), k(NULL), WaveNo(NULL), 
		    BandWidth(NULL), istart(NULL), iend(NULL),
		    iCO(NULL), iCO2(NULL), iH2O(NULL), iMix(NULL),
		    Tn(NULL), k_CO(NULL), k_CO2(NULL), k_H2O(NULL), 
		    k2_CO(NULL), k2_CO2(NULL), k2_H2O(NULL), dT(ZERO),
		    Ninterp(0), MixType(SNBCK_OVERLAP_OPTICALLY_THIN),
		    EvalType(SNBCK_EVAL_ONLINE)
{}




/*********************************************************************
 * Destructor                                                        *
 *********************************************************************/
SNBCK :: ~SNBCK() {  Deallocate(); }


/*********************************************************************
 * Array allocators/deallocators                                     *
 *********************************************************************/
void SNBCK :: AllocateQuad() {

  // deallocate just in case
  DeallocateQuad();

  // loop over bands
  if (Nquad>0) {
    g = new double [Nquad];
    w = new double [Nquad];
  }
}

void SNBCK :: AllocateBands() {

  // deallocate just in case
  DeallocateBands();

  // intialize
  if (Nbands>0) {
    WaveNo    = new double[Nbands];
    BandWidth = new double[Nbands];
    istart = new int[Nbands];
    iend   = new int[Nbands];
    iCO    = new int[Nbands];
    iCO2   = new int[Nbands];
    iH2O   = new int[Nbands];
    iMix   = new int[Nbands];
  }

}

void SNBCK :: AllocateInterp() {

  // deallocate just in case
  DeallocateInterp();

  // intialize (Nquad should be uniform - uncorrelated case)
  if (Nbands>0 && Nquad>0 && Ninterp>0) {
    Tn     = new double[Ninterp];
    k_CO   = new double**[Nbands];
    k_CO2  = new double**[Nbands];
    k_H2O  = new double**[Nbands];
    k2_CO  = new double**[Nbands];
    k2_CO2 = new double**[Nbands];
    k2_H2O = new double**[Nbands];
    for (int i=0; i<Nbands; i++) {
      k_CO[i]   = new double*[Nquad];
      k_CO2[i]  = new double*[Nquad];
      k_H2O[i]  = new double*[Nquad];
      k2_CO[i]  = new double*[Nquad];
      k2_CO2[i] = new double*[Nquad];
      k2_H2O[i] = new double*[Nquad];
      for (int j=0; j<Nquad; j++) {
	k_CO[i][j]   = new double[Ninterp];
	k_CO2[i][j]  = new double[Ninterp];
	k_H2O[i][j]  = new double[Ninterp];
	k2_CO[i][j]  = new double[Ninterp];
	k2_CO2[i][j] = new double[Ninterp];
	k2_H2O[i][j] = new double[Ninterp];
      } /* endfor */
    } /* endfor */
  } /* endif */

}

void SNBCK :: AllocateAbs() {

  // deallocate just in case
  DeallocateAbs();

  if (Nquad>0 && Nbands>0) {
    k = new double*[Nbands];
    for (int i=0; i<Nbands; i++) {
      k[i] = new double [Nquad];
    }/* endfor */
  }/* endif */

}


void SNBCK :: DeallocateQuad() {

  if (g != NULL) delete[] g; g=NULL;
  if (w != NULL) delete[] w; w=NULL;
}



void SNBCK :: DeallocateBands() {
  if (WaveNo    != NULL)  delete[] WaveNo;       WaveNo = NULL;
  if (BandWidth != NULL)  delete[] BandWidth; BandWidth = NULL;
  if (istart    != NULL)  delete[] istart;       istart = NULL;
  if (iend      != NULL)  delete[] iend;           iend = NULL;
  if (iCO       != NULL)  delete[] iCO;             iCO = NULL;
  if (iCO2      != NULL)  delete[] iCO2;           iCO2 = NULL;
  if (iH2O      != NULL)  delete[] iH2O;           iH2O = NULL;
  if (iMix      != NULL)  delete[] iMix;           iMix = NULL;
}


void SNBCK :: DeallocateInterp() {
  if (k_CO != NULL) { 

    for (int i=0; i<Nbands; i++) {
      for (int j=0; j<Nquad; j++) {
	delete[] k_CO  [i][j];
	delete[] k_CO2 [i][j];
	delete[] k_H2O [i][j];
	delete[] k2_CO [i][j];
	delete[] k2_CO2[i][j];
	delete[] k2_H2O[i][j];
      } /* endfor */
      delete[] k_CO  [i];
      delete[] k_CO2 [i];
      delete[] k_H2O [i];
      delete[] k2_CO [i];
      delete[] k2_CO2[i];
      delete[] k2_H2O[i];
    } /* endfor */
    delete[]   k_CO;  k_CO  =NULL;
    delete[]  k_CO2;  k_CO2 =NULL;
    delete[]  k_H2O;  k_H2O =NULL;
    delete[]  k2_CO;  k2_CO =NULL;
    delete[] k2_CO2;  k2_CO2=NULL;
    delete[] k2_H2O;  k2_H2O=NULL;

  }/* endif */

  if (Tn != NULL)  delete[] Tn;  Tn = NULL;
}


void SNBCK :: DeallocateAbs() {

  if (k != NULL) { 
    for (int i=0; i<Nbands; i++) delete[] k[i];
    delete[] k; k=NULL;
  }
}


void SNBCK :: Deallocate() {  
  DeallocateQuad();    DeallocateBands(); 
  DeallocateInterp();  DeallocateAbs();
}


/*********************************************************************
 * SNBCK :: Setup                                                    *
 *                                                                   *
 * Load the dataset EM2C of Soufiani and Taine (1997) for CO, CO2,   *
 * and H2O which is used to compute the SNB model parameters.        *
 * Setup the narrow/wide bands, determine quadrature rule, and       *
 * allocate absorbsion coefficient array.  For precalculated case,   *
 * generate spline interpolants for reference absorbsion coefficient.*
 *********************************************************************/
//------------------------------------------------
// for precalculation of cummulative distribution 
// function 
//------------------------------------------------
void SNBCK :: Setup( const int quad_type,  // quadrature type
		     const int quad_points,// number of quadrature points
		     const int Nlump,      // number of narrow bands lumped together
		     const bool optimize,  // attempt to optimize band lumping (Nlump>1)
		     const char *PATH,     // Current path
		     const double p_ref,    // pressure [atm]
		     const double xco_ref,  // mole fraction oh CO
		     const double xh2o_ref, // mole fraction oh H2O
		     const double xco2_ref, // mole fraction oh CO2
		     const double xo2_ref,  // mole fraction of O2
		     const int Nint )       // number of interpolation points 

{
  // deallocate just to be sure
  Deallocate();

  // store eval flag
  EvalType = SNBCK_EVAL_PRECALC;

  // only uncorrelated mixture rule can be used
  MixType = SNBCK_OVERLAP_UNCORRELATED;

  // get SNB model parameters for gas mixture (units: cm, atm, K)
  SNB.LoadParams(PATH);

  // lump bands together
  SetupBands( Nlump, optimize );

  // setup quadrature
  SetupQuad( quad_type , quad_points );

  // we are precalculating the absorbsion coefficient
  PreCalculateAbsorb( p_ref, xco_ref, xh2o_ref, xco2_ref, xo2_ref, Nint);

  // allocate absorbsion coefficient
  AllocateAbs();

}


//------------------------------------------------
// for online inversion of cummulative 
// distribution function 
//------------------------------------------------
void SNBCK :: Setup( const int quad_type,  // quadrature type
		     const int quad_points,// number of quadrature points
		     const int Nlump,      // number of narrow bands lumped together
		     const bool optimize,  // attempt to optimize band lumping (Nlump>1)
		     const char *PATH,     // Current path
		     const int mix_rule )  // mixture rule

{

  // deallocate just to be sure
  Deallocate();

  // store eval flag
  EvalType = SNBCK_EVAL_ONLINE;

  // store mixture rule
  MixType = mix_rule;

  // get SNB model parameters for gas mixture (units: cm, atm, K)
  SNB.LoadParams(PATH);

  // lump bands together
  SetupBands( Nlump, optimize );

  // setup quadrature
  SetupQuad( quad_type, quad_points );

  // allocate absorbsion coefficient
  AllocateAbs();

}


/*********************************************************************
 * SetupQuad                                                         *
 *                                                                   *
 * Given the quadrature type, setup the weights and values at each   *
 * point. Gauss-Legendre (Gauss to Dr. Groth) is used.               *
 *********************************************************************/
void SNBCK :: SetupQuad( const int quad_type,    // quadrature type
			 const int quad_points ) // number of quadrature points
{
  // declares
  int cnt;

  // allocate necessary arrays
  Nquad = quad_points;
  AllocateQuad();

  // get quadrature (g defined on 0 to 1)
  switch (quad_type) {
  case GAUSS_LEGENDRE:
    gauleg(0, 1, g, w, quad_points);
    break;
  default:
    cerr << "SNBCK::SetupQuad(): Invalid value for quadrature type." << endl;
    exit(-1);
  } /* end switch */


  
//------------- Use half of Gauss points defined on [-1,1] -------------//
//   // allocate temporary storage
//   double* xx = new double[2*quad_points];
//   double* ww = new double[2*quad_points];

//   // get quadrature (g defined on 0 to 1)
//   switch (quad_type) {
//   case GAUSS_LEGENDRE:
//     gauleg(-1, 1, xx, ww, 2*quad_points);
//     break;
//   default:
//     cerr << "SNBCK::SetupQuad(): Invalid value for quadrature type." << endl;
//     exit(-1);
//   } /* end switch */

//   // copy arrays
//   for (int i=0; i<quad_points; i++){
//     g[i] = xx[i+quad_points];
//     w[i] = ww[i+quad_points];
//   }

//   // clean up
//   delete[] xx; delete[] ww;
//----------------------------------------------------------------------//

}


/*********************************************************************
 * LumpBands                                                         *
 *                                                                   *
 * Lump uniform narrow SNB bands into using the uniformly spaced or  *
 * optimized band lumping strategy.  The latter is the selective     *
 * regrouping of narrow bands based on the radiative properties of   *
 * CO2. Basically, transparent bands of CO2 are lumped into wider    *
 * bands than emitting CO2 bands.  Note that there are only two      *
 * optimized strategies implementd and that they are essentially     *
 * 'Hard-coded'.  See                                                *
 *       Goutiere, Charette, Kiss, Numerical Heat Transfer - Part B, *
 *       v41, pp. 361-381. (2002)                                    *
 * and                                                               *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, AIAA-99-3679 (1999)     *
 *********************************************************************/
void SNBCK :: SetupBands( const int Nlump,      // number of narrow bands lumped together
			  const bool optimize ) // attempt to optimize band lumping (Nlump>1)
{ 
	       
  //------------------------------------------------
  // Optimized band lumping, Liu et al. (1999)
  //------------------------------------------------
  if ( optimize && Nlump>=30 ) {

    // Optimized regrouping based on the radiative properties of CO2. 
    // This gas possesses four large primary bands separated by semitransparent bands. 
    Nbands = 9;
    AllocateBands();

    // set wide band spectral indexes
    istart[ 0] =   0;  iend[ 0] =  11; // 12 transparent
    istart[ 1] =  12;  iend[ 1] =  42; // 31 emitting
    istart[ 2] =  43;  iend[ 2] =  71; // 29 transparent
    istart[ 3] =  72;  iend[ 3] =  92; // 21 emitting
    istart[ 4] =  93;  iend[ 4] = 125; // 33 transparent
    istart[ 5] = 126;  iend[ 5] = 146; // 21 emitting
    istart[ 6] = 147;  iend[ 6] = 181; // 35 transparent
    istart[ 7] = 182;  iend[ 7] = 204; // 23 emitting
    istart[ 8] = 205;  iend[ 8] = 366; // 162 transparent


  //------------------------------------------------
  // Optimized band lumping, Goutiere et al. (2002)
  //------------------------------------------------
  } else if ( optimize && Nlump<30 && Nlump>1 ) {

    // Optimized regrouping based on the radiative properties of CO2. 
    // This gas possesses four large primary bands separated by semitransparent bands. 
    // On the emission bands, seven or eight narrow bands can be regrouped, whereas
    // the regrouping on the more transparent bands can be of 12 to 40.
    Nbands = 21;
    AllocateBands();

    // set wide band spectral indexes
    istart[ 0] =   0;  iend[ 0] =  11; // 12 transparent
    istart[ 1] =  12;  iend[ 1] =  19; // 8 emitting
    istart[ 2] =  20;  iend[ 2] =  27; // 8 emitting
    istart[ 3] =  28;  iend[ 3] =  35; // 8 emitting
    istart[ 4] =  36;  iend[ 4] =  42; // 7 emitting
    istart[ 5] =  43;  iend[ 5] =  71; // 29 transparent
    istart[ 6] =  72;  iend[ 6] =  78; // 7 emitting
    istart[ 7] =  79;  iend[ 7] =  85; // 7 emitting
    istart[ 8] =  86;  iend[ 8] =  92; // 7 emitting
    istart[ 9] =  93;  iend[ 9] = 125; // 33 transparent
    istart[10] = 126;  iend[10] = 132; // 7 emitting
    istart[11] = 133;  iend[11] = 139; // 7 emitting
    istart[12] = 140;  iend[12] = 146; // 7 emitting
    istart[13] = 147;  iend[13] = 181; // 35 transparent
    istart[14] = 182;  iend[14] = 189; // 8 emitting
    istart[15] = 190;  iend[15] = 197; // 8 emitting
    istart[16] = 198;  iend[16] = 204; // 7 emitting
    istart[17] = 205;  iend[17] = 244; // 40 transparent
    istart[18] = 245;  iend[18] = 284; // 40 transparent
    istart[19] = 285;  iend[19] = 324; // 40 transparent
    istart[20] = 325;  iend[20] = 366; // 42 transparent


  //------------------------------------------------
  // Regular uniformly lumped bands, Liu et al. (1999)
  //------------------------------------------------
  } else {

    // compute the number of new bands
    // Note that the remaining bands at the higher frequencies are dropped
    // when there are narrow bands left over.
    Nbands = int(SNB.Nbands/Nlump);
    AllocateBands();
    
    // loop over the lumped bands
    istart[0] = 0;
    for (int i=0; i<Nbands-1; i++) {
      iend[i]     = istart[i]+Nlump-1;
      istart[i+1] = iend[i]+1;
    } /* endfor */
    iend[Nbands-1] = istart[Nbands-1]+Nlump-1;

  } /* endif */ 

  //------------------------------------------------
  // End strategy specific
  //------------------------------------------------


  // Determine if any wide bands overlap.
  // iMix also holds the total number of overlapping 
  // gases for a particular band.
  for (int i=0; i<Nbands; i++) {
    iH2O[i] = 0; iCO2[i] = 0; iCO[i] = 0; iMix[i] = 0;
    for (int n=istart[i]; n<=iend[i]; n++) {
      iH2O[i] += SNB.liH2O[n];
      iCO2[i] += SNB.liCO2[n];
      iCO[i]  += SNB.liCO[n];
      iMix[i] += ( iH2O[i]>0 ? 1 : 0 ) + 
	         ( iCO2[i]>0 ? 1 : 0 ) + 
	         ( iCO [i]>0 ? 1 : 0 );
    } /* endfor */
  } /* endfor - original bands */



  // determine the new bandwidth centers
  double vmin;
  for (int i=0; i<Nbands; i++) {

    // compute the left bracket of the band
    vmin = SNB.WaveNo[istart[i]] - SNB.BandWidth[istart[i]]/TWO;

    // determine the new bandwidth
    BandWidth[i] = ZERO;
    for (int n=istart[i]; n<=iend[i]; n++) BandWidth[i] += SNB.BandWidth[n];

    // compute the new bandwidth center
    WaveNo[i] = vmin + BandWidth[i]/TWO;

  } /* endfor */

}


/*********************************************************************
 * SNBCK :: CalculateAbsrob                                          *
 *                                                                   *
 * Wrapper function.  Calculates absorbsion coefficient directly     *
 * or using spline interpolants built from the reference state.      *
 *********************************************************************/
void SNBCK :: CalculateAbsorb( const double p,        // pressure [atm]
			       const double T,        // temperature [K]
			       const double xco,      // mole fraction oh CO
			       const double xh2o,     // mole fraction oh H2O
			       const double xco2,     // mole fraction oh CO2
			       const double xo2,      // mole fraction of O2
			       const double xsoot )   // volume fraction of soot  
{
  
  switch (EvalType) {

    // if we are doing online inversion
  case SNBCK_EVAL_ONLINE:
    CalculateAbsorb_Direct( p, T, xco, xh2o, xco2, xo2, xsoot ); 
    break;

    // if we are precalculating the absorbsion coefficient
  case SNBCK_EVAL_PRECALC:
    CalculateAbsorb_Interp( p, T, xco, xh2o, xco2, xo2, xsoot ); 
    break;

  default:
    cerr << "Error in SNBCK::CalculateAbsorb() : Invalid "
	 << "value for evaluation type flag.\n";
    exit(-1);

  }

}



/*********************************************************************
 * SNBCK :: CalculateAbsrob_Direct                                   *
 *                                                                   *
 * Given the temperature, pressure, species composition (CO, CO2,    *
 * H2O, soot), quadrature rule, the absorbsion coefficient for each  *
 * band, at each quadrature point, are computed.  The dataset EM2C   *
 * of Soufiani and Taine (1997) for CO, CO2, and H2O is used to      *
 * compute the SNB model parameters.  For overlapping bands, the     *
 * optically thin approximation proposed by Liu et al. (2001) is     * 
 * used.  Bands are lumped together using the strategy proposed by   *
 * Liu et al. (1999).                                                *
 *                                                                   *
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/
void SNBCK :: CalculateAbsorb_Direct( const double p,        // pressure [atm]
				      const double T,        // temperature [K]
				      const double xco,      // mole fraction oh CO
				      const double xh2o,     // mole fraction oh H2O
				      const double xco2,     // mole fraction oh CO2
				      const double xo2,      // mole fraction of O2
				      const double xsoot )   // volume fraction of soot  
{

  // get SNB model parameters for gas mixture (units: cm, atm, K)
  SNB.ComputeSNB( p, T, xco, xh2o, xco2, xo2 );

 
  //------------------------------------------------
  // Compute absorbsion coefficient for each wide  
  // band at each quadrature point
  //------------------------------------------------

  switch (MixType) {


  //-----------------------------------------------
  // OPTICALLY THIN APPROXIMATION proposed by Liu et al. (2001)
  //-----------------------------------------------
  case SNBCK_OVERLAP_OPTICALLY_THIN:
    
    // loop over each wide band
    for (int v=0; v<Nbands; v++) {
      
      // loop over quadrature points
      for (int i=0; i<Nquad; i++) {

	// if an absorbing gray band, then compute the absorbsion coeffient,
	// else its zero
	if ( (iH2O[v] && xh2o>MICRO) || 
	     (iCO2[v] && xco2>MICRO) ||
	     (iCO[v] && xco>MICRO) )
	  k[v][i] = AbsorptionCoeffSNBCK( g[i], SNB.B_Thin, SNB.S_Thin, 
					  SNB.liMix, istart[v], iend[v] );
	else k[v][i]=ZERO;

      } /* end for - quadrature */
    } /* end for - band */


    break;

  //-----------------------------------------------
  // OPTICALLY THICK APPROXIMATION proposed by Lacis and Oinas. (1991)
  //-----------------------------------------------
  case SNBCK_OVERLAP_OPTICALLY_THICK:
    
    // loop over each wide band
    for (int v=0; v<Nbands; v++) {
      
      // loop over quadrature points
      for (int i=0; i<Nquad; i++) {

	// if an absorbing band, then compute the absorbsion coeffient,
	// else its zero
	if ( (iH2O[v] && xh2o>MICRO) || 
	     (iCO2[v] && xco2>MICRO) ||
	     (iCO[v] && xco>MICRO) )
	  k[v][i] = AbsorptionCoeffSNBCK( g[i], SNB.B_Thick, SNB.S_Thick, 
					  SNB.liMix, istart[v], iend[v] );
	else k[v][i]=ZERO;

      } /* end for - quadrature */
    } /* end for - band */


    break;


  //-----------------------------------------------
  // UNCORRELATED APPROXIMATION proposed by Liu et al. (2001)
  // Fully CORRELATED see Liu et al. (2001)
  //-----------------------------------------------
  case SNBCK_OVERLAP_UNCORRELATED:

    // loop over each wide band
    for (int v=0; v<Nbands; v++) {
      

      // loop over quadrature points 
      for (int i=0; i<Nquad; i++) {

	// add each active component
	k[v][i] = ZERO;
	if (iH2O[v] && xh2o>MICRO) 
	  k[v][i] += AbsorptionCoeffSNBCK( g[i], SNB.B_H2O, SNB.S_H2O, 
					   SNB.liH2O, istart[v], iend[v] );
	if (iCO2[v] && xco2>MICRO) 
	  k[v][i] += AbsorptionCoeffSNBCK( g[i], SNB.B_CO2, SNB.S_CO2, 
					   SNB.liCO2, istart[v], iend[v] );
	if (iCO[v] && xco>MICRO) 
	  k[v][i] += AbsorptionCoeffSNBCK( g[i],  SNB.B_CO,  SNB.S_CO,  
					   SNB.liCO, istart[v], iend[v] );

      } /* end for - quadrature */
    } /* end for - bands */
    break;
        

  //-----------------------------------------------
  // shouldn't get here
  //-----------------------------------------------
  default:
    cerr << "SNBCK::SetupSNB: Invalid value for mix_rule." << endl;
    exit(-1);
    break;
    
  } /* end switch */


  //-----------------------------------------------
  // End model specific
  //-----------------------------------------------


  // Add soot component to absorbsion coefficient [cm^-1]
  for (int v=0; v<Nbands; v++) {
    for (int i=0; i<Nquad; i++) {
      k[v][i] += 5.5 * xsoot * WaveNo[v];
    } /* end for */
  } /* end for */


}


/*********************************************************************
 * PreCalculateAbsorb                                                *
 *                                                                   *
 * Given the reference state (pressure, 1 atm, and species comp.     *
 * (CO, CO2, H2O, soot) ), and the quadrature rule,                  *
 * the cummulative distribution function for each quadrature point   *
 * is precalculated (curvefit) wrt temperature. This procedure       *
 * follows that proposed by Liu and Smallwood (2004).  This requires *
 * the assumption that the SNB model parameter B is independant of   *
 * species concentration.  The absorbsion coefficients for each      *
 * species computed from the curvefits are then treated as           *
 * uncorrelated (i.e. added together).  The dataset EM2C             *
 * of Soufiani and Taine (1997) for CO, CO2, and H2O is used to      *
 * compute the SNB model parameters.  For overlapping bands, the     *
 * optically thin approximation proposed by Liu et al. (2001) is     * 
 * used.  Bands are lumped together using the strategy proposed by   *
 * Liu et al. (1999).                                                *
 *                                                                   *
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/
void SNBCK :: PreCalculateAbsorb( const double p_ref,    // pressure [atm]
				  const double xco_ref,  // mole fraction oh CO
				  const double xh2o_ref, // mole fraction oh H2O
				  const double xco2_ref, // mole fraction oh CO2
				  const double xo2_ref,  // mole fraction of O2
				  const int Nint )       // number of interpolation points 
{

  //------------------------------------------------
  // Initialize
  //------------------------------------------------

  // make sure enough interpolation points
  if (Nint<3) Ninterp = SNB.Npoints;
  else Ninterp = Nint;

  // allocate interpolation parameters
  AllocateInterp();


  // for uniform spacing
  double Tmin = SNB.Tmin;
  double Tmax = SNB.Tmax;
  dT = (Tmax-Tmin)/double(Ninterp-1);

  //------------------------------------------------
  // loop over each temperature point
  //------------------------------------------------
  for (int n=0; n<Ninterp; n++) {

    // compute temperautre
    Tn[n] = Tmin + double(n)*dT;

    // get SNB model parameters for gas mixture (units: cm, atm, K)
    // and normalize by the reference state
    SNB.ComputeRefSNB( p_ref, Tn[n], xco_ref, xh2o_ref, xco2_ref, xo2_ref );

      
    //------------------------------------------------
    // loop over each wide band, inverting the 
    // cummulative distribution function
    //------------------------------------------------
    for (int v=0; v<Nbands; v++) {
      for (int i=0; i<Nquad; i++) {
	   
	// add each active component
	k_CO[v][i][n] = ZERO; k_CO2[v][i][n] = ZERO; k_H2O[v][i][n] = ZERO;
	if (iH2O[v])
	  k_H2O[v][i][n] += AbsorptionCoeffSNBCK( g[i], SNB.B_H2O, SNB.S_H2O, 
						  SNB.liH2O, istart[v], iend[v] );
	if (iCO2[v])
	  k_CO2[v][i][n] += AbsorptionCoeffSNBCK( g[i], SNB.B_CO2, SNB.S_CO2, 
						  SNB.liCO2, istart[v], iend[v] );
	if (iCO[v])
	  k_CO[v][i][n] += AbsorptionCoeffSNBCK( g[i],  SNB.B_CO,  SNB.S_CO,  
						 SNB.liCO, istart[v], iend[v] );
      } 
    } 
    //------------------------------------------------
    // endfor - bands, quad points
    //------------------------------------------------


  }
  //------------------------------------------------
  // endfor - Temperature
  //------------------------------------------------


  // Build a natural cubic spline
  for (int v=0; v<Nbands; v++) 
    for (int i=0; i<Nquad; i++) {
      spline( Tn,  k_CO[v][i], Ninterp, 1.1E30, 1.1E30,  k2_CO[v][i] );
      spline( Tn, k_CO2[v][i], Ninterp, 1.1E30, 1.1E30, k2_CO2[v][i] );
      spline( Tn, k_H2O[v][i], Ninterp, 1.1E30, 1.1E30, k2_H2O[v][i] );
    }

}

/*********************************************************************
 * SNBCK :: CalculateSNB_Interp                                      *
 *                                                                   *
 * Given the temperature, pressure, species composition (CO, CO2,    *
 * H2O, soot), quadrature rule, the absorbsion coefficient for each  *
 * band, at each quadrature point, are computed.  The dataset EM2C   *
 * of Soufiani and Taine (1997) for CO, CO2, and H2O is used to      *
 * compute the SNB model parameters.  Here, the values are computed  *
 * using a predertimened cubic spline. This procedure                *
 * follows that proposed by Liu and Smallwood (2004).  This requires *
 * the assumption that the SNB model parameter B is independant of   *
 * species concentration.  The absorbsion coefficients for each      *
 * species computed from the curvefits are then treated as           *
 * uncorrelated (i.e. added together).  The dataset EM2C             *
 * of Soufiani and Taine (1997) for CO, CO2, and H2O is used to      *
 * compute the SNB model parameters.  For overlapping bands, the     *
 * optically thin approximation proposed by Liu et al. (2001) is     * 
 * used.  Bands are lumped together using the strategy proposed by   *
 * Liu et al. (1999).                                                 *                                                                   *
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/
void SNBCK :: CalculateAbsorb_Interp( const double p,        // pressure [atm]
				      const double T,        // temperature [K]
				      const double xco,      // mole fraction oh CO
				      const double xh2o,     // mole fraction oh H2O
				      const double xco2,     // mole fraction oh CO2
				      const double xo2,      // mole fraction of O2
				      const double xsoot )   // volume fraction of soot  
{

  // check to make sure 
  if (k_CO == NULL) {
    cerr << "Error in SNBCK::CalculateAbsorb_Interp() : Need to run "
	 << "SNBCK::PreCalculateAbsorb() first.\n";
    exit(-1);
  }

  // declares
  double kk_CO, kk_CO2, kk_H2O;

  for (int v=0; v<Nbands; v++) 
    for (int i=0; i<Nquad; i++) {

      // interpolate using the cubic spline
      splint( Tn,  k_CO[v][i],  k2_CO[v][i], dT, Ninterp, T,  kk_CO );
      splint( Tn, k_CO2[v][i], k2_CO2[v][i], dT, Ninterp, T, kk_CO2 );
      splint( Tn, k_H2O[v][i], k2_H2O[v][i], dT, Ninterp, T, kk_H2O );

      // compute absorbsion coefficient using uncorrelated approximation
      // See Liu and Smallwood (2004)
      k[v][i] = kk_CO*xco + kk_CO2*xco2 + kk_H2O*xh2o;
      k[v][i] *= p;
      
      // add soot
      k[v][i] += 5.5 * xsoot * WaveNo[v];

    } /* endfo */


}


/*********************************************************************
 * SNBCK :: BandAverage                                              *
 *                                                                   *
 * Computes the band avereaged value where phi is an array of length *
 * Nquad containing the values of phi at each quadrature point.      *
 *********************************************************************/
double SNBCK :: BandAverage( const double *phi, const int v )  {
  double avg = ZERO;
  for (int i=0; i<Nquad; i++) avg += w[i]*phi[i];
  return avg;
}


/*********************************************************************
 * SNBCK :: Transmissivity                                           *
 *                                                                   *
 * Computes the band avereaged gas transmissivity for an isothermal  *
 * and homogeneous path.  See Liu et al. (2001).                     *
 *      L - path length in cm                                        *
 *********************************************************************/
double SNBCK :: Transmissivity( const double L, const int v ){

  double tau=ZERO;
  for (int i=0; i<Nquad; i++) tau += w[i]*exp(-k[v][i]*L);
  return tau;

}





/*********************************************************************
 ******************************** VARIOUS FUNCTIONS ******************
 *********************************************************************/


/*********************************************************************
 * g, g_prime                                                        *
 *                                                                   *
 * Computes the cummulative distribution function given the          *
 * absorbsion coefficient, k, and the SNB model parameters, B, S.    *
 * Also, the derivative w.r.t. k is computed.                        *
 * See:                                                              *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 *********************************************************************/
double g( const double k, const double B, const double S ){
  // check for optically thick case
  //   if ( B > GRAY_GAS_CUTOFF)  {
  //     if (k>=S) return ONE;
  //     else return ZERO;
  //   } /* endif */


  // compute g1 using Eq. (7) in Liu et al. (2000)
  double a = HALF*sqrt( PI*B*S ); 
  double b = HALF*sqrt( PI*B/S ); 
  double temp = HALF * erfc( a/sqrt(k) - b*sqrt(k) ) +
                HALF * erfc( a/sqrt(k) + b*sqrt(k) )*exp(PI*B);
  return temp;
}

double g_prime( const double k, const double B, const double S ){

  // check for optically thick case
  //   if ( B > GRAY_GAS_CUTOFF )  return MILLION;

  // compute dg/dk using Eq. (7) in Liu et al. (2000)
  double a = HALF*sqrt( PI*B*S );
  double b = HALF*sqrt( PI*B/S );
  double c1 = -pow( a/sqrt(k)-b*sqrt(k), 2 );
  double c2 = -pow( a/sqrt(k)+b*sqrt(k), 2 );
  double d1 = -HALF*a*pow(k,-THREE/TWO) - HALF*b*pow(k,-HALF);
  double d2 = -HALF*a*pow(k,-THREE/TWO) + HALF*b*pow(k,-HALF);
  double temp = -HALF*( exp(c1)*d1 + exp(c2+PI*B)*d2 );
  temp *= TWO/sqrt(PI);
  return temp;
}

void g_lumped( const double*B, const double*S,  const int*iFlag,
	       const int Nstart, const int Nend, const int Nlump,
	       const double k, double &gg, double &dgg){

  
  gg = ZERO;  dgg = ZERO; 
  for (int n=Nstart; n<=Nend; n++) {
    if (iFlag[n]) {
      // compute g1 using Eq. (7) in Liu et al. (2000)
      gg += g(k,B[n],S[n]) / Nlump;
      // the derivative of g
      dgg += g_prime(k,B[n],S[n])/Nlump;
    }
  } /* endfor */

}

/*********************************************************************
 * getAbsorptionCoeff                                                *
 *                                                                   *
 * Computes the absorbsion coefficient based on the SNB model        *
 * parameters input.  Newton-Raphson iteration combined with         *
 * bisectionis used (See Numerical Recipes in C, 2nd ed, 1992).      *
 * As an initial guess, the value of k where the laplace function    *
 * f peaks is used.  This procedure is outlined in                   *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 *********************************************************************/
double AbsorptionCoeffSNBCK( const double g_val, const double*B, 
			     const double*S,  const int*iFlag,
			     const int Nstart, const int Nend){
  
  // set iteration paramters
  int Nmax = 200;    // max iterations
  double TOL = 1E-6; // error tolerance
  
  // the function values
  double gn;                 // g(k)
  double dgn;                // g'(k)
  double kn,k0;
  double dkn, dk0; 
  double kmin = ZERO;        // left bracket
  double kmax = MILLION;     // right bracket
  double gmin = ZERO-g_val;  // function value at left bracket
  double gmax = ONE-g_val;   // function value at right bracket

  // number of bands to lump together
  int Nlump = (Nend - Nstart) + 1;


  // check to make sure that we are not on the bounds
  if (fabs(gmin)==ZERO) return kmin;
  if (fabs(gmax)==ZERO) return kmax;


  //------------------------------------------------
  // Initial Guess
  //------------------------------------------------

  // determine the value of k where f(k) peaks for initial guess
  // depending upon how many bands are lumped together, take an average
  double k_peak = ZERO;
  for (int n=Nstart; n<=Nend; n++) {
    if (iFlag[n]) k_peak += PeakAbsorptionCoeffSNBCK(B[n], S[n])/Nlump;
  }

  // make sure guess is within bracket to start
  if (k_peak>kmin && k_peak<kmax) kn = k_peak;
  else kn = (kmin+kmax)/TWO;

  // depending upon how many bands are lumped together,
  // compute the cummulative dist function and its derivative
  g_lumped( B, S, iFlag, Nstart, Nend, Nlump, kn, gn, dgn);
  gn -= g_val; 

  // initialize the previous stepsizes
  dkn = dk0 = fabs(kmax - kmin);

  // search is already oriented such that g(kmin)<0


  //------------------------------------------------
  // Newton-Raphson iteration
  //------------------------------------------------

  int i; //counter
  for ( i=1; i<=Nmax; i++){

    // use bisection if Newton out of range or not decreasing 
    // fast enough
    if ( (((kn-kmax)*dgn-gn)*((kn-kmin)*dgn-gn) > ZERO) || 
	 (fabs(TWO*gn) > fabs(dk0*dgn)) ) {
      dk0 = dkn;
      dkn = HALF*(kmax-kmin);
      kn = kmin + dkn;

    // Newton acceptable 
    } else {
      dk0 = dkn;
      dkn = gn/dgn;
      kn -= dkn;
    }


    // Convergence test
    if (fabs(dkn)/kn<TOL || dkn!=dkn) break;      
      
    // depending upon how many bands are lumped together,
    // compute the cummulative dist function and its derivative
    g_lumped( B, S, iFlag, Nstart, Nend, Nlump, kn, gn, dgn);
    gn -= g_val; 

    // change bisection bracket
    if ( gn < ZERO)  kmin=kn;
    else kmax = kn;
   

  } /* end Newton-Raphson */
  


  //------------------------------------------------
  // Convergence
  //------------------------------------------------

  // check to see if it converged, if not, print an error
  if ( i>=Nmax || dkn!=dkn ) {
    cerr << "getAbsorptionCoeff(): Newton-Raphson iteration failed to converge." << endl;
    cerr << "getAbsorptionCoeff(): err = " << fabs(dkn) 
	 << ", f(x) = " << gn << ", f'(x) = " << dgn << ", x = " << g_val     
	 << ", kmin = " << kmin << ", k = " << kn << ", kmax = " << kmax <<  endl;
    cerr << "getAbsorptionCoeff(): k_peak = " << k_peak << ", it = " << i <<  endl;
    cerr << "B = " << B[Nstart] << ", S = " << S[Nstart] <<  endl;
    cerr << endl;
  }


  // return the value
  return kn;

}


/*********************************************************************
 * getPeakAbsorptionCoeff                                            *
 *                                                                   *
 * Computes the peak absorbsion coefficient where the value of the   *
 * laplace transform f, Eq (5) in Liu et al. (2000), is a maximum.   *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 * Assumes the gas follows a Malkmus distribution.                   *
 *********************************************************************/
double PeakAbsorptionCoeffSNBCK( const double B, const double S ){

  // NOTE:
  // -> See Eq. (19) in Lacias and Oinas, J Geophysical Research, v96, 
  // pp. 9027-9063, 1991.
  double k1 = sqrt( pow(PI*B/THREE,2) + ONE ) - ONE; 
  k1 *= THREE*S / (PI*B);

  // -> return peak, since b>0, then k2<0 always  
  if ( k1<ZERO || k1!=k1 ) {
    cerr << "getPeakAbsorptionCoeff(): Error computing k_peak." << endl;
    cerr << "getPeakAbsorptionCoeff(): B = " << setw(18) << B 
	 << ", S = " << setw(18) << S 
      	 << ", k_peak = " << k1
	 <<  endl;
  }
  return k1;


}


/*********************************************************************
 * printDistFunction                                                 *
 *                                                                   *
 * Outputs the cummulative distribution function, g, Eq (7) in:      *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 * Assumes the gas follows a Malkmus distribution.                   *
 *********************************************************************/
void OutputDistFunction(  const double B, const double S, 
			  const double kmin, const double kmax,
			  const char *file ){

  // declares
  int Nk = int(1E6);
  double delta_k = (kmax-kmin)/double(Nk-1);
  int Ng = 100;
  double delta_g = ONE/double(Ng-1);
  double k1, g1, Ba[1]={B}, Sa[1]={S};
  int iFlag[1];


  // open output stream and print the header
  ofstream out(file,ios::out);
  out << "VARIABLES = \"k [cm^-1]\" \"g\""
      << endl;

  //------------------------------------------------
  // Evaluating g(k)
  //------------------------------------------------
  out << "ZONE T=\"Evaluating g(k)\"" << endl;

  // compute
  for ( int i=0; i<Nk; i++){

    // compute g1 using Eq. (7) in Liu et al. (2000)
    k1 = kmin + i*delta_k;
    g1 = g(k1,B,S);

    // output
    out << setw(24) << scientific << setprecision(15) << k1 
	<< setw(24) << scientific << setprecision(15) << g1 << endl;
  }  

  //------------------------------------------------
  // Evaluating k(g)
  //------------------------------------------------
  out << endl << "ZONE T=\"Evaluating k(g)\"" << endl;

  // compute
  for ( int i=0; i<Ng; i++){

    // compute k1 inverting Eq. (7) in Liu et al. (2000)
    g1 = i*delta_g;
    k1 = AbsorptionCoeffSNBCK( g1, Ba, Sa, iFlag, 0, 0);

    // output
    out << setw(24) << scientific << setprecision(15) << k1 
	<< setw(24) << scientific << setprecision(15) << g1 << endl;
  }  


  // close output stream
  out.close();

}



/*********************************************************************
 * getTransmissivity                                                 *
 *                                                                   *
 * Computes the transmissivity over an isothermal and homogeneous    *
 * path assuming a Malkmus gas. Eq (4) in:                           *
 *      F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass       *
 *      Trans., v43, 2000. pp. 3119-3135.                            *
 *                                                                   *
 *      L - path length in cm                                        *
 *********************************************************************/
double TransmissivitySNB( const double B, const double S, const double L ){

  double Bnew = B;
  if (Bnew<PICO) Bnew = PICO;
    double temp1 = sqrt( ONE + FOUR*S*L/(PI*B) ) - ONE;
  double temp2 = exp( -PI*B/TWO * temp1 );
  return temp2;

}


/*********************************************************************
 * LineOfSightIntens                                                 *
 *                                                                   *
 * Compute the ratio I2/I1 where I1 is the knwn, initial intensity at*
 * a point x=0 and I2 is the downstream intensity along a line of    *
 * sight at x=L
 *                                                                   *
 *      L - path length in cm                                        *
 *********************************************************************/
double LineOfSightIntens( const double I1, const double Ib, 
			  const double k, const double L ){

  return  ( I1 * exp(-k*L) + 
	    Ib*( ONE - exp(-k*L)) );

}

double LineOfSightIntens( const double I1, const double Ib, 
			  const double tau ){

  return  ( I1*tau + Ib*(ONE-tau) );

}


/*********************************************************************
 * Convert fortran "D" double notation to cpp "E" notation.  Returns *
 * a double.                                                         * 
 *********************************************************************/
double FortranToCppD( string str ) {
  str.replace(str.find("D"),1,"e");
  return atof(str.c_str());
}


/********************************************************
 * Approximation of error function                      *
 ********************************************************/
double erfc_a(const double x) {
  static double p =   0.3275911;
  static double a1 =  0.254829592;
  static double a2 = -0.284496736;
  static double a3 =  1.421413741;
  static double a4 = -1.453152027;
  static double a5 =  1.061405429;
  double t = ONE / (ONE+p*x);
  return ( a1*t + a2*pow(t,2) + a3*pow(t,3) + a4*pow(t,4) + a5*pow(t,5) ) *
          exp( -pow(x,2) );
}

double erfc_exp(const double x1, const double x2) {
  static double p =   0.3275911;
  static double a1 =  0.254829592;
  static double a2 = -0.284496736;
  static double a3 =  1.421413741;
  static double a4 = -1.453152027;
  static double a5 =  1.061405429;
  double t = ONE / (ONE+p*x1);
  return ( a1*t + a2*pow(t,2) + a3*pow(t,3) + a4*pow(t,4) + a5*pow(t,5) ) *
          exp( -pow(x1,2) + x2 );
}




// /*********************************************************************
//  * functions to return numerical quadrature data                     *
//  *                                                                   *
//  * Given the quadrature type, setup return weights and values at each*
//  * point. Quadrature taken from :                                    *
//  *   F. Liu, G.J. Smallwood, O.L. Gulder, Int. J. Heat Mass Trans.,  *
//  *   v43, 2000. pp. 3119-3135.                                       *
//  *********************************************************************/
// /* Total number of quad points */
// int nquad( const int quad_type ) {

//   switch (quad_type){

//     // 1 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_1P:
//       return 1;
    
//     // 2 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_2P:
//       return 2;

//     // 4 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_4P:
//       return 4;

//     // 4 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_4P:
//       return 4;

//     // 5 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_5P:
//       return 5;

//     // 7 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_7P:
//       return 7;

//     // default case, shouldn't get here
//     default:
//       cerr << "Nquad(): Invalid quadrature type." << endl;
//       exit(-1);
//       break;
//   }

// }

// /* return value and weight at index i */
// void quad( const int quad_type, const int i, double &g, double &w  ) {

//   switch (quad_type){

//     // 1 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_1P:
//       switch (i){
//       case (0):
//       default:
// 	g = 0.57735;  w = 1.00000;
// 	break;
//       } /* endswitch */
//       break;

//     // 2 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_2P:
//       switch (i){
//       case (0):
// 	g = 0.33998;  w = 0.65215;
// 	break;
//       case (1):
//       default:
// 	g = 0.86114;  w = 0.34785;
// 	break;
//       } /* endswitch */
//       break;

//     // 4 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_4P:
//       switch (i){
//       case (0):
// 	g = 0.18343;  w = 0.36268;
// 	break;
//       case (1):
// 	g = 0.52553;  w = 0.31371;
// 	break;
//       case (2):
// 	g = 0.79667;  w = 0.22238;
// 	break;
//       case (3):
//       default:
// 	g = 0.96029;  w = 0.10123;
// 	break;
//       } /* endswitch */
//       break;

//     // 4 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_4P:
//       switch (i){
//       case (0):
// 	g = 0.20930;  w = 0.41246;
// 	break;
//       case (1):
// 	g = 0.59170;  w = 0.34112;
// 	break;
//       case (2):
// 	g = 0.87174;  w = 0.21070;
// 	break;
//       case (3):
//       default:
// 	g = 1.00000;   w = 0.03571;
// 	break;
//       } /* endswitch */
//       break;

//     // 5 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_5P:
//       switch (i){
//       case (0):
// 	g = 0.16528;  w = 0.32754;
// 	break;
//       case (1):
// 	g = 0.47792;  w = 0.29204;
// 	break;
//       case (2):
// 	g = 0.73877;  w = 0.22489;
// 	break;
//       case (3):
// 	g = 0.91953;  w = 0.13331;
// 	break;
//       case (4):
//       default:
// 	g = 1.00000;  w = 0.02222;
// 	break;
//       } /* endswitch */
//       break;
      
//     // 7 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_7P:
//       switch (i){
//       case (0):
// 	g = 0.00000;  w = 0.04500;
// 	break;
//       case (1):
// 	g = 0.15541;  w = 0.24500;
// 	break;
//       case (2):
// 	g = 0.45000;  w = 0.32000;
// 	break;
//       case (3):
// 	g = 0.74459;  w = 0.24500;
// 	break;
//       case (4):
// 	g = 0.90000;  w = 0.05611;
// 	break;
//       case (5):
// 	g = 0.93551;  w = 0.05125;
// 	break;
//       case (6):
//       default:
// 	g = 0.98449;  w = 0.03764;
// 	break;
//       } /* endswitch */
//       break;

//     // default case, shouldn't get here
//     default:
//       cerr << "Quad(): Invalid quadrature type." << endl;
//       exit(-1);

//   } /* endswitch */

// }

// /* return number of points, and the value and weight arrays */
// void quad( const int quad_type, int &n, double* (&g), double* (&w)  ) {

//   switch (quad_type){

//     // 1 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_1P:
//       n = 1;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.57735;
//       w[0] = 1.00000;
//       break;
    
//     // 2 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_2P:
//       n = 2;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.33998;  w[0] = 0.65215;
//       g[1] = 0.86114;  w[1] = 0.34785;
//       break;

//     // 4 Point Gauss-Legendre quadrature scheme
//     case GAUSS_LEGENDRE_4P:
//       n = 4;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.18343;  w[0] = 0.36268;
//       g[1] = 0.52553;  w[1] = 0.31371;
//       g[2] = 0.79667;  w[2] = 0.22238;
//       g[3] = 0.96029;  w[3] = 0.10123;
//       break;

//     // 4 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_4P:
//       n = 4;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.20930;  w[0] = 0.41246;
//       g[1] = 0.59170;  w[1] = 0.34112;
//       g[2] = 0.87174;  w[2] = 0.21070;
//       g[3] = 1.00000;  w[3] = 0.03571;
//       break;

//     // 5 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_5P:
//       n = 5;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.16528;  w[0] = 0.32754;
//       g[1] = 0.47792;  w[1] = 0.29204;
//       g[2] = 0.73877;  w[2] = 0.22489;
//       g[3] = 0.91953;  w[3] = 0.13331;
//       g[4] = 1.00000;  w[4] = 0.02222;
//       break;

//     // 7 Point Gauss-Lobatto quadrature scheme
//     case GAUSS_LOBATTO_7P:
//       n = 7;
//       g = new double[n];
//       w = new double[n];
//       g[0] = 0.00000;  w[0] = 0.04500;
//       g[1] = 0.15541;  w[1] = 0.24500;
//       g[2] = 0.45000;  w[2] = 0.32000;
//       g[3] = 0.74459;  w[3] = 0.24500;
//       g[4] = 0.90000;  w[4] = 0.05611;
//       g[5] = 0.93551;  w[5] = 0.05125;
//       g[6] = 0.98449;  w[6] = 0.03764;
//       break;

//     // default case, shouldn't get here
//     default:
//       cerr << "SNBCK::SetupQuadrature: Invalid quadrature type." << endl;
//       exit(-1);
//       break;
//   }

// }

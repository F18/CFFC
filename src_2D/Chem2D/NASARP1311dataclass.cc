 /********************* NASARP1311class.cc **********************************
   This class pasrses the data files and contains constructors to extract
   useful thermodynamic and transport data as a function of temperature.

  NOTES:  The string compare function works differently 
          depending if its complied with cxx or g++????
          hence the use of -D_GNU_GCC_296 compiler flag.

        assosiated header file: 
                  NASARP1311.h

*****************************************************************************/

#ifndef _NASARP1311_DATA_INCLUDED
#include "NASARP1311data.h"
#endif // _NASARP1311_DATA_INCLUDED

/**************************************************************************
********************* THERMODATA CLASS CONSTRUCTORS ***********************
***************************************************************************/

/********************* GET DATA FOR SPECIES *******************************
  This constructor reads in all the thermodynamic and transport data for
  a species "spec" which is a string for the species, ie. O, O2, CH4, etc.

  It is rather obfuscated as it is parsing FORTRAN formated data files and 
  has to read everything in as strings, then convert them to doubles and 
  ints depending on there spacings as outlined in the NASA RP1311 Report
  in the appendices, pg73 for thermo.inp, pg96 for trans.inp.
***************************************************************************/

//NASARP1311data::NASARP1311data(string &spec)
void NASARP1311data::Getdata(string spec, const char *PATH, const int &trans_data)
{

  // set transport flag
  trans_type = trans_data;

  //Set Data path for thermo.inp and trans.inp
  Set_Path_Names(PATH);

  // get the thermodynamic data
  GetThermoData(spec);

  // get the tranport data and initialize function pointers
  if (trans_type == TRANSPORT_NASA) {
    pt_Viscosity = &NASARP1311data::Viscosity_NASA;
    pt_ThermalConduct = &NASARP1311data::ThermalConduct_NASA;
    pt_dViscositydT = &NASARP1311data::dViscositydT_NASA;
    GetTransDataNASA(spec);
  } else if (trans_type == TRANSPORT_LENNARD_JONES) {
    pt_Viscosity = &NASARP1311data::Viscosity_LJ_Poly;            // fit a polynomial
    pt_ThermalConduct = &NASARP1311data::ThermalConduct_LJ_Poly;  // fit a polynomial
    //pt_Viscosity = &NASARP1311data::Viscosity_LJ;
    //pt_ThermalConduct = &NASARP1311data::ThermalConduct_LJ;
    pt_dViscositydT = &NASARP1311data::dViscositydT_LJ;
    GetTransDataLJ(spec);   
  } else {
    cerr << "Error - Getdata(): Never heard of this transport dataset!";
    exit(-1);
  }
   
}

/********************* GET THERMO DATA *******************************/
void NASARP1311data::GetThermoData(const string spec) {
  
  if(temp_ranges_thermo != 0) { deallocate_thermo_data(); }


  //------------ READ IN WHAT SPECIES ---------------------//
  species = spec;
  //add a space at the end of the string, necessary for search accuracy
  species += ' ';

  //variables
  string type;
  int spec_l;
  
  spec_l = species.length();

  //------------- OPEN DATA FILES ------------------------//
  ifstream thermodatafile(datafilename_thermo);

  //Check to see if successful
  if(thermodatafile.fail()){ 
    cerr<<"\nError opening file: "<<datafilename_thermo<<endl;
    exit(0); 
  }	 

  if(thermodatafile.fail()){ 
    cerr<<"\nError opening file: "<<datafilename_trans<<endl;
    exit(0); 
  }

  //------------------------------------------------------------//
  //----------- FIND SPECIES THERMODYNAMIC DATA ----------------//
  //------------------------------------------------------------//
  // sort through file until desired species
  do{
    getline(thermodatafile,type);
  } 
#ifdef _GNU_GCC_296 //redhat 7.3 gcc-2.96 
  while( type.compare(species,0,spec_l) != 0 && !thermodatafile.eof() );
#else //alpha & redhat 9.0 gcc-3.22
  while( type.compare(0,spec_l,species) != 0 && !thermodatafile.eof() );
#endif

  if(thermodatafile.eof()){
    cerr<<" Species: "<<species<<" not in "<<datafilename_thermo<<endl;
    exit(1);
  }

  //---------- GET SPECIES DATA -------------------------------//
  // number of temperature intervals  
  getline(thermodatafile,type);

  string T = type.substr(0,2); 
  temp_ranges_thermo = stoi(T);
  

  //get Molecular Weight and Heat of Formation
  string M = type.substr(53,12);
  mol_mass = stof(M)/1000.0; //convert from g/mol to kg/mol
  string H = type.substr(68,12);
  heatofform = stof(H)/mol_mass; //convert from J/mol to J/kg

  //container for cofficients  
  thermo_data = new thermocoef[temp_ranges_thermo]; 

  double Cp_coef[9]; //1-7 +b1,b2
  
  //---------- Loop through temperature intervals -------------//     
  //get each range of coefficient data
  for( int i = 0; i<temp_ranges_thermo; i++){

    //line 3
    getline(thermodatafile,type);
    //string Temp = type.substr(2,12);
    string Temp = type.substr(1,10);
    thermo_data[i].Low_range_in(stof(Temp));

    //Temp = type.substr(13,21);
    Temp = type.substr(13,9);
    thermo_data[i].High_range_in(stof(Temp));

    Temp = type.substr(22,1);
    //thermo_data[i].Numexponents_in(stoi(Temp)); 

    //convert from J/mol to J/kg by / my mol_mass
    Temp = type.substr(66,80);
    thermo_data[i].DeltaHref_in(stof(Temp)/mol_mass);

    // Everything is set up around 7 exponents so check!
    //if(thermo_data[i].Numexponents() != 7){
    //  cerr<<thermo_data[i].Numexponents()<<" exponents for Cp, should be 7";
    //}

    //line4
    getline(thermodatafile,type);
  
    //Get 7 polynomial coefficients, messy but it works.
    for(int j=0; j<5; j++){
      //convert Fortran "D" notation to C "e" notation
      //Temp = type.substr(j*16,16+j*16);
      Temp = type.substr(j*16,16);
      string::size_type pos = Temp.find("D");
      Temp.replace(pos,1,"e");
      Cp_coef[j] = stof(Temp); 

      //cout.precision(12);
      //cout<<Cp_coef[j]<<endl; 
    }
  
    //line5
    getline(thermodatafile,type);

    for(int j=0; j<5; j++){
      //convert Fortran "D" notation to C "e" notation
      if( j < 2){
	//Temp = type.substr(j*16,16+j*16);
	Temp = type.substr(j*16,16);
	string::size_type pos = Temp.find("D");
	Temp.replace(pos,1,"e");
	Cp_coef[j+5] = stof(Temp);

	//cout.precision(12);
	//cout<<Cp_coef[j+5]<<endl;
      }
      //skips the 3rd entry in the 5th line as it
      //is not used. ie. its 0.0000D+00
      else if (j > 2){
	//Temp = type.substr(j*16,16+j*16);
	Temp = type.substr(j*16,16);
	string::size_type pos = Temp.find("D");
	Temp.replace(pos,1,"e");
	Cp_coef[j+4] = stof(Temp);
	
	//cout.precision(12);
	//cout<<Cp_coef[j+4]<<endl;
      }
    }
    
    //move into thermocoef class
    thermo_data[i].Thermo_coef_in(Cp_coef);
  }

  //close thermodynamic data file
  thermodatafile.close();

}

/********************* GET TRANS DATA *******************************/
void NASARP1311data::GetTransDataNASA(const string spec) {

  if(temp_ranges_V != 0) { deallocate_trans_V(); }
  if(temp_ranges_C != 0) { deallocate_trans_C(); }
  
  //------------ READ IN WHAT SPECIES ---------------------//
  species = spec;
  //add a space at the end of the string, necessary for search accuracy
  species += ' ';
  
  //variables
  string line;
  int spec_l;
   
  spec_l = species.length();
  
  
  //------------- OPEN DATA FILES ------------------------//
  ifstream transdatafile(datafilename_trans);
  
  //Check to see if successful
  if(transdatafile.fail()){ 
    cerr<<"\nError opening file: "<<datafilename_trans<<endl;
    exit(0); 
  }

  //------------------------------------------------------------//
  //-------------- FIND SPECIES TRANSPORT DATA -----------------//
  //------------------------------------------------------------//

  //adding 20 spaces to make sure the data is for pure species
  //not binary data which would have another species name 
  //in spaces 17-31. 
  species += "                    "; //20 spaces
  spec_l = species.length();
 
  //check each line until proper species, or end of file.
  do{
    getline(transdatafile,line);
  } 
#ifdef _GNU_GCC_296 //redhat 7.3 gcc-2.96 
  while( line.compare(species,0,spec_l) != 0 && !transdatafile.eof() );
#else //alpha & redhat 9.0 gcc-3.22
  while( line.compare(0,spec_l,species) != 0 && !transdatafile.eof() );
#endif

  if(transdatafile.eof()){
    cerr<<" Species: "<<species<<" not in "<<datafilename_trans<<endl;
    exit(1);
  }
  
  //---------- GET SPECIES DATA -------------------------------//
  // number of temperature intervals for viscosity (V) and
  // thermal conductivity (C)  
  string V,C;

  if( line.substr(34,1) == "V"){
    V = line.substr(35,1);
    temp_ranges_V = stoi(V);
  }
  if( line.substr(36,1) == "C"){
    C = line.substr(37,1);
    temp_ranges_C = stoi(C);
  }

  //container for coefficients of size temp_ranges
   trans_viscosity = new transcoef[temp_ranges_V]; 
   trans_thermconduct = new transcoef[temp_ranges_C];

  double V_coef[4]; //A,B,C,D see pg 96 NASA RP1311 Vol. 2 
  double C_coef[4];

  //---------- Loop through Viscosity values -------------//     
  //get each range of coefficient data
  for( int i = 0; i<temp_ranges_V; i++){

    //line 2,3
    getline(transdatafile,line);
    
    //Viscosity (V)
    if(line.substr(1,1) == "V"){ 
      
      //lower and upper temp ranges for data
      //string Temp = line.substr(3,8);
      string Temp = line.substr(2,7);
      trans_viscosity[i].Low_range_in(stof(Temp));

      //Temp = line.substr(11,8);
      Temp = line.substr(11,7);
      trans_viscosity[i].High_range_in(stof(Temp));
 
      //data values
      for(int j=0; j<4; j++){
	Temp = line.substr(20+j*15,15);

	//fix missing E+ in data files
	if( Temp.find("E ") > 0 && Temp.find("E ") < 15){
	  Temp.replace(Temp.find("E "),2,"E+");
	}
	V_coef[j] = stof(Temp);
       
	//cout.precision(12);
	//cout<<V_coef[j]<<endl;
      } 
      //move into transocoef class
      trans_viscosity[i].Trans_coef_in(V_coef); 
    }
  }
 
  //--------- Loop through Thermal Conductivity values ---------//
  //get each range of coefficient data
  for( int i = 0; i<temp_ranges_C; i++){
    
    //line 4,5
    getline(transdatafile,line);
    
    //Thermal Conductivity (V)
    if(line.substr(1,1) == "C"){ 
      
      //lower and upper temp ranges for data
      //string Temp = line.substr(3,8);
      string Temp = line.substr(2,7);
      trans_thermconduct[i].Low_range_in(stof(Temp));

      //Temp = line.substr(11,8);
      Temp = line.substr(11,7);
      trans_thermconduct[i].High_range_in(stof(Temp));
   
      //data values
      for(int j=0; j<4; j++){
	Temp = line.substr(20+j*15,15);	

	//fix missing E+ in data files
	if( Temp.find("E ") > 0 && Temp.find("E ") < 15){
	  Temp.replace(Temp.find("E "),2,"E+");
	}
	C_coef[j] = stof(Temp);
       
	//cout.precision(12);
	//cout<<C_coef[j]<<endl;
      } 
      //move into transocoef class
      trans_thermconduct[i].Trans_coef_in(C_coef); 
    }
  }
  //fix up species name
  species = species.replace(species.find("                     "),21,""); //21 spaces
  //close transport data file
  transdatafile.close();
}

/********************* GET TRANS DATA *******************************/
void NASARP1311data::GetTransDataLJ(const string spec) {
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
  in.open(datafilename_trans);
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
  
  // get the potential well-depth / boltzman constant
  eps = stof( line.substr(20,10) ) * BOLTZMANN;
 
  // get the collision diameter in Angstroms and convert to m
  sigma = stof( line.substr(30,10) ) * 1.0E-10;  
  
  // get the dipole moment in Debye and convert to m^1.5 J^0.5
  mu = stof( line.substr(40,10) ) * 1.0E-25*sqrt(10.0); 
  
  // get the polarizability in Angstroms^3 and convert to m^3
  alpha = stof( line.substr(50,10) ) * 1.0E-30; 
 
  // get the rotational relaxation parameter
  Zrot = stof( line.substr(60,10) );
  
  // close the transport data file
  in.close();

  //cout << "\n Species "<< species <<": g = " << g << "\t eps = " << eps
  //     << "\t sigma = " << sigma << "\n\t mu = " << mu  
  //     << "\t alpha = " << alpha << "\t Zrot = " << Zrot;;


  //----------------------------------
  //    Polynomial fitting
  //----------------------------------
  ofstream outfile;

  // polynomial parameters
  const int np = 100,   // number of interpolation points
            degree = 3; // degree of polynomial

  // counter for term
  int ndeg = 0;        

  // temperature interval parameters
  double dT, T,          // stepsize, temperature point
         Tmin = 200.0,   // min value
         Tmax = 3500.0;  // max value

  // temporary arrays
  double *w, *Tlog, *spvisc, *spcond;

  // assume atmospheric pressure   ?????
  double P = 101325.0;

  // deallocate, just to be sure
  if(temp_ranges_V != 0) { deallocate_trans_V(); }
  if(temp_ranges_C != 0) { deallocate_trans_C(); }

  //container for coefficients of size temp_ranges (one range)
  trans_viscosity = new transcoef[1]; 
  trans_thermconduct = new transcoef[1];
  double V_coef[degree+1];  
  double C_coef[degree+1];
  
  // allocate temporary storage
  w = new double[np];
  Tlog = new double[np];
  spvisc = new double[np];
  spcond = new double[np];

  //double Tlog[np], spvisc[np], spcond[np], w[np];

  //
  // generate array of log(T) values
  //
  dT = (Tmax - Tmin)/double(np-1);     
  for (int n = 0; n < np; ++n) {
    T = Tmin + dT*double(n);
    Tlog[n] = log(T);

    spvisc[n] = log( Viscosity_LJ(T) );
    spcond[n] = log( ThermalConduct_LJ(T, P) );
    w[n] = -1.0;
  }
  
  // call the polynomial fitting function
  polyfit(np, Tlog, spvisc, w, degree, ndeg, 0.0, V_coef);
  polyfit(np, Tlog, spcond, w, degree, ndeg, 0.0, C_coef);

  // move into transcoef class
  trans_viscosity[0].Trans_coef_in(V_coef);
  trans_thermconduct[0].Trans_coef_in(C_coef);

  //
  // evaluate max fit errors for viscosity
  //
  double val, fit, err, relerr;
  outfile.open("Viscosity.dat", ios::out);
  for (int n = 0; n < np; ++n) {
    val = exp(spvisc[n]);
    fit = exp( poly3(Tlog[n], V_coef) ); 
    err = fit - val;
    relerr = err/val;
    outfile << exp(Tlog[n]) << "  " << val << "  " << fit 
	    << "  " << relerr << endl;
  }
  outfile.close();
 
  //
  // evaluate max fit errors for conductivity
  //
  outfile.open("Thermal_Conductivity.dat", ios::out);
  for (int n = 0; n < np; ++n) {
    val = exp(spcond[n]);
    fit = exp( poly3(Tlog[n], C_coef) ); 
    err = fit - val;
    relerr = err/val;
    outfile << exp(Tlog[n]) << "  " << val << "  " << fit 
	    << "  " << relerr << endl;
  }
  outfile.close();

  // deallocate memory
  delete[]  w;
  delete[]  Tlog;
  delete[]  spvisc;
  delete[]  spcond;
 
}


/**************** Find coefficient to use *************************************/
//for use with thermodynamic data
int NASARP1311data::Which_coef(double &Temp){

  //check temperature and find proper coef[i] range to use, find which i 
  //find which set of coefficients corresponds to temp { 0... temp_ranges }
  for(int i=0; i<temp_ranges_thermo; i++){
    if (Temp >= thermo_data[i].Low_range() && Temp <= thermo_data[i].High_range()){   
      return i; 
    }
  }
  
  //check if below range and if it is use lowest data set.
  if (Temp < thermo_data[0].Low_range() && Temp > 0.0 ){
    Temp = thermo_data[0].Low_range();
    return 0;   
  }
  
  //check if too high for range, temperature is negative which just doesn't work, or its a "nan"
  else if (Temp > thermo_data[temp_ranges_thermo -1].High_range() || Temp < 0.0 || Temp != Temp ){
    cerr<<"\n Requested Temperature of "<<Temp<<" for "<<species<<" out of range"<<endl;
    exit(1);
  }


  /****** ORIGINAL ******/
//   //check if below range and if it is use lowest data set.
//   if (Temp < thermo_data[0].Low_range() && Temp > 0.0 ){
//     Temp = thermo_data[0].Low_range();
//     return 0;   
//   }  
//   //check if too high for range, temperature is negative which just doesn't work, or its a "nan"
//   else if (Temp > thermo_data[temp_ranges_thermo -1].High_range() || Temp < 0.0 || Temp != Temp ){
//     cerr<<"\n Requested Temperature of "<<Temp<<" for "<<species<<" out of range"<<endl;
//     exit(1);
//   }
//   //so it must be in the range of the data set thus 
//   //find which set of coefficients corresponds to temp { 0... temp_ranges }
//   else{
//     for(int i=0; i<=temp_ranges_thermo; i++){
//       if (Temp >= thermo_data[i].Low_range() && Temp <= thermo_data[i].High_range()){    
// 	return i;
//       }
//     }
//   }

}

/**************** Find coefficient to use *************************************/
// for use with transport data
int NASARP1311data::Which_coef(double &Temp,
			       transcoef *trans,
			       int temp_ranges){

  //check temperature and find proper coef[i] range to use, find which i 
  //find which set of coefficients corresponds to temp { 0... temp_ranges }  
  for(int i=0; i<temp_ranges; i++){
    if (Temp >= trans[i].Low_range() && Temp <= trans[i].High_range()){  
      return i;
    }
  }

  //check if below range and if it is use lowest data set.
  if (Temp < trans[0].Low_range() && Temp > 0.0 ){
    Temp = trans[0].Low_range();
    return 0;   
  }
  
  //check if to high for range or temperature is negative which just doesn't work in Kelvin
  else if (Temp > trans[temp_ranges -1].High_range() || Temp < 0.0){
    cerr<<"\n Requested Temperature of "<<Temp<<" for "<<species<<" out of range"<<endl;
    exit(1);
  }

//   /****** ORIGINAL ******/
//   //check if below range and if it is use lowest data set.
//   if (Temp < trans[0].Low_range() && Temp > 0.0 ){
//     Temp = trans[0].Low_range();
//     return 0; 
//   }
//   //check if to high for range or temperature is negative which just doesn't work in Kelvin
//   else if (Temp > trans[temp_ranges -1].High_range() || Temp < 0.0){
//     cerr<<"\n Requested Temperature of "<<Temp<<" for "<<species<<" out of range"<<endl;
//     exit(1);
//   }
//   //so it must be in the range of the data set thus 
//   //find which set of coefficients corresponds to temp { 0... temp_ranges }
//   else{
//     for(int i=0; i<=temp_ranges; i++){
//       if (Temp >= trans[i].Low_range() && Temp <= trans[i].High_range()){   
// 	return i;
//       }
//     }
//   }

}


/********************** Heat Capacity ******************************************/
double NASARP1311data::HeatCapacity_p(double Temp){
 
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp);
  
  //use polynomial and return Cp. Ref. NASA RP1311 Appendix A
//   double Cp = R_UNIVERSAL *  ( thermo_data[i].Thermo_coef(0)/(Temp*Temp)
// 		     + thermo_data[i].Thermo_coef(1)/(Temp)
// 		     + thermo_data[i].Thermo_coef(2)
// 		     + thermo_data[i].Thermo_coef(3)*Temp
// 		     + thermo_data[i].Thermo_coef(4)*Temp*Temp
// 		     + thermo_data[i].Thermo_coef(5)*Temp*Temp*Temp
// 		     + thermo_data[i].Thermo_coef(6)*Temp*Temp*Temp*Temp);

  // form requiring less computation 
  double Cp = R_UNIVERSAL * ( (thermo_data[i].Thermo_coef(0)/Temp
		     + thermo_data[i].Thermo_coef(1))/Temp
		     + thermo_data[i].Thermo_coef(2)
		     + Temp*(thermo_data[i].Thermo_coef(3)
		     + Temp*(thermo_data[i].Thermo_coef(4)			     
		     + Temp*(thermo_data[i].Thermo_coef(5)
		     + Temp*thermo_data[i].Thermo_coef(6)))));


  return Cp/mol_mass;  
}


/********************** Enthalpy  **********************************************/
double NASARP1311data::Enthalpy(double Temp){
  
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp);

  //use polynomial and return Cp. Ref. NASA RP1311 Appendix A
//   double H = R_UNIVERSAL*Temp * ( - thermo_data[i].Thermo_coef(0)/(Temp*Temp)
// 			+ log(Temp)*thermo_data[i].Thermo_coef(1)/(Temp)
// 			+ thermo_data[i].Thermo_coef(2)
// 			+ thermo_data[i].Thermo_coef(3)*Temp/2.0
// 			+ thermo_data[i].Thermo_coef(4)*Temp*Temp/3.0
// 			+ thermo_data[i].Thermo_coef(5)*Temp*Temp*Temp/4.0
// 			+ thermo_data[i].Thermo_coef(6)*Temp*Temp*Temp*Temp/5.0 
// 			+ thermo_data[i].Thermo_coef(7)/Temp );

  // form requiring less computation 
  double H = R_UNIVERSAL * ( - thermo_data[i].Thermo_coef(0)/Temp
			     + log(Temp)*thermo_data[i].Thermo_coef(1)
			     + Temp*(thermo_data[i].Thermo_coef(2)
			     + Temp*(thermo_data[i].Thermo_coef(3)/TWO
			     + Temp*(thermo_data[i].Thermo_coef(4)/THREE	     
			     + Temp*(thermo_data[i].Thermo_coef(5)/FOUR
			     + Temp*(thermo_data[i].Thermo_coef(6)/FIVE))))) 
 			     + thermo_data[i].Thermo_coef(7) );	     
  // J/kg
  // h = H + DeltaHref - heatofform
  //cout<<"\nH "<<H/mol_mass<<" "<<DeltaHref()<<" "<<heatofform<<" "<<(H/mol_mass + DeltaHref() - heatofform)<<endl; 
  return (H/mol_mass + DeltaHref() - heatofform);  
 
}

/********************** Enthalpy  **********************************************/
double NASARP1311data::Enthalpy_mol(double Temp){
  
  int i = Which_coef(Temp);
  //use polynomial and return Cp. Ref. NASA RP1311 Appendix A
  double H = R_UNIVERSAL*Temp * ( - thermo_data[i].Thermo_coef(0)/(Temp*Temp)
			+ log(Temp)*thermo_data[i].Thermo_coef(1)/(Temp)
			+ thermo_data[i].Thermo_coef(2)
			+ thermo_data[i].Thermo_coef(3)*Temp/2.0
			+ thermo_data[i].Thermo_coef(4)*Temp*Temp/3.0
			+ thermo_data[i].Thermo_coef(5)*Temp*Temp*Temp/4.0
			+ thermo_data[i].Thermo_coef(6)*Temp*Temp*Temp*Temp/5.0 
			+ thermo_data[i].Thermo_coef(7)/Temp );
  // J/mol
  return (H); //  + (DeltaHref() - heatofform)*mol_mass);
}

/********************** Enthalpy Derivate (dH/dT) *******************************/
double NASARP1311data::Enthalpy_prime(double Temp){
  return HeatCapacity_p(Temp);  
}

/********************** Entropy ***********************************************/
double NASARP1311data::Entropy(double Temp){
 
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp);
 
  //use polynomial and return Cp. Ref. NASA RP1311 Appendix A
  double S = R_UNIVERSAL * ( - thermo_data[i].Thermo_coef(0)/(Temp*Temp*2.0)
		   - thermo_data[i].Thermo_coef(1)/(Temp)
		   + thermo_data[i].Thermo_coef(2)*log(Temp)
		   + thermo_data[i].Thermo_coef(3)*Temp
		   + thermo_data[i].Thermo_coef(4)*Temp*Temp/2.0
		   + thermo_data[i].Thermo_coef(5)*Temp*Temp*Temp/3.0
		   + thermo_data[i].Thermo_coef(6)*Temp*Temp*Temp*Temp/4.0 
		   + thermo_data[i].Thermo_coef(8) );
  //J/kgK
  return S/mol_mass;  
}

/********************** Entropy ***********************************************/
double NASARP1311data::Entropy_mol(double Temp){
  //J/molK
  return Entropy(Temp)*mol_mass;
}



/********************** Viscosity *********************************************/
// compute the viscosity using NASA polynomial
double NASARP1311data::Viscosity_NASA(double Temp){
  
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp,trans_viscosity,temp_ranges_V);

  //use polynomial and return Cp. Ref. NASA RP1311 Appendix E
  double eta =   trans_viscosity[i].Trans_coef(0)*log(Temp)
               + trans_viscosity[i].Trans_coef(1)/Temp
               + trans_viscosity[i].Trans_coef(2)/(Temp*Temp)
               + trans_viscosity[i].Trans_coef(3);              
 
  //converted from mircoposie (dynes*s/cm^2)x10^-6 to kg/(m*s)
  // 1 (dynes*s/cm^2) = .1 kg/(m*s)
  return exp(eta)*1e-7;
}

//The derivative of viscosity w.r.t Temperature using NASA polynomial
double NASARP1311data::dViscositydT_NASA(double Temp){
  
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp,trans_viscosity,temp_ranges_V);
 
  //use polynomial and return Cp. Ref. NASA RP1311 Appendix E
  double eta =   trans_viscosity[i].Trans_coef(0)*log(Temp)
               + trans_viscosity[i].Trans_coef(1)/Temp
               + trans_viscosity[i].Trans_coef(2)/(Temp*Temp)
               + trans_viscosity[i].Trans_coef(3);              
 
  double prod =   trans_viscosity[i].Trans_coef(0)/Temp
                  - trans_viscosity[i].Trans_coef(1)/(Temp*Temp)
                 -2.0*trans_viscosity[i].Trans_coef(2)/(Temp*Temp*Temp);
            
  //converted from mircoposie (dynes*s/cm^2)x10^-6 to kg/(m*s)
  // 1 (dynes*s/cm^2) = .1 kg/(m*s)
  return exp(eta)*1e-7*prod;
}

/****************** Thermal Conductivity ***************************************/
double NASARP1311data::ThermalConduct_NASA(double Temp, double Press){
  
  //find appropriate coefficients for temperature range
  int i = Which_coef(Temp,trans_thermconduct,temp_ranges_C);

  //use polynomial and return Cp. Ref. NASA RP1311 Appendix E
  double eta =   trans_thermconduct[i].Trans_coef(0)*log(Temp)
               + trans_thermconduct[i].Trans_coef(1)/Temp
               + trans_thermconduct[i].Trans_coef(2)/(Temp*Temp)
               + trans_thermconduct[i].Trans_coef(3);              
  
  //converted from (W/cm*K)x10^-6 to W/m*K
  return exp(eta)*1e-4;
}

// /******************* Plot out curve fits **************************************
//   Use Gnuplot to plot out the polynomial representations of the data for
//   Cp, H, S, mu, and k depening on the handle passed in.  Uses the gnuplot
//   ternary operator to plot each of the temperature ranges on the same figure,
//   much like a piecewise function.  

//   This constructor makes use of stringstreams which are very useful for 
//   constructing strings out of many types of data.  They are used in the 
//   same way as cout and fstreams are used. ie. the "<<" and ">>" operators.

//   The output is a postscript file named after the species and the physical 
//   property plotted.
// ******************************************************************************/
// void NASARP1311data::Plot_data(string type){
  
//   //number of ranges
//   int num_range; // could be temp_ranges_thermo, temp_ranges_C, temp_ranges_V
//   string units;  // units for each type

//   // Create gnuplot script file 
//   const char gnu_script[80] = "plot.temp";
//   ofstream gnu_out(gnu_script);
 
//   //Check to see if successful
//   if(gnu_out.fail()){ 
//     cerr<<"\nError opening file: "<<gnu_script<<endl;
//     exit(1); 
//   }	   

//   //string streams for ease of setting up equations and ranges
//   ostringstream *equation = new ostringstream[temp_ranges_thermo]; 
//   ostringstream *segment = new ostringstream[temp_ranges_thermo];
//   ostringstream range;

//   //-----------------------------------------------------------------//
//   //-------- Set up polynomial equation -----------------------------//
//   //-----------------------------------------------------------------//
 
//   //---------------- Heat Capacity ---------------------------------//
//   if(type.compare("Cp") == 0){

//     num_range = temp_ranges_thermo;
//     //plot all temperature ranges on one plot
//     for(int i=0; i<num_range; i++){

//       //specific polynomial
//       equation[i]<<R_UNIVERSAL/mol_mass<<"*("<<thermo_data[i].Thermo_coef(0)<<"/(x*x) +"
// 		 <<thermo_data[i].Thermo_coef(1)<<"/x +"
// 		 <<thermo_data[i].Thermo_coef(2)<<"+"
// 		 <<thermo_data[i].Thermo_coef(3)<<"*x +"
// 		 <<thermo_data[i].Thermo_coef(4)<<"*x**2 +"
// 		 <<thermo_data[i].Thermo_coef(5)<<"*x**3 +"
// 		 <<thermo_data[i].Thermo_coef(6)<<"*x**4 )"	
// 		 <<" : 1/0 ";
      
//       //piecewise range using ternary opperator
//       segment[i]<<thermo_data[i].Low_range()<<"<=x && x<="<<thermo_data[i].High_range()<<" ? "; 
//     }
//     //overall range 
//     range<<" ["<<0<<":"<<thermo_data[num_range-1].High_range()<<"] ";
//     units = "(J/(kg*K))";
//   }

//   //------------------ Enthalpy -------------------------------------//
//   else if(type.compare("H") == 0){

//     num_range = temp_ranges_thermo;
//     //plot all temperature ranges on one plot
//     for(int i=0; i<num_range; i++){

//       //specific polynomial
//       equation[i]<<R_UNIVERSAL/mol_mass<<"*x*("<<-thermo_data[i].Thermo_coef(0)<<"/(x*x) +"
// 		 <<thermo_data[i].Thermo_coef(1)<<"*log(x)/x +"
// 		 <<thermo_data[i].Thermo_coef(2)<<"+"
// 		 <<thermo_data[i].Thermo_coef(3)<<"*x/2.0 +"
// 		 <<thermo_data[i].Thermo_coef(4)<<"*x**2/3.0 +"
// 		 <<thermo_data[i].Thermo_coef(5)<<"*x**3/4.0 +"
// 		 <<thermo_data[i].Thermo_coef(6)<<"*x**4/5.0 +"
// 	         <<thermo_data[i].Thermo_coef(7)<<"/x ) "
// 		 <<" : 1/0 ";
      
//       //piecewise range using ternary opperator 
//       segment[i]<<thermo_data[i].Low_range()<<"<=x && x<="<<thermo_data[i].High_range()<<" ? "; 
//     }
//     //overall range 
//     range<<" ["<<0<<":"<<thermo_data[num_range-1].High_range()<<"] ";  
//     units = "(J/kg)"; 
//   }

//   //------------------ Entropy -------------------------------------//
//   else if(type.compare("S") == 0){
    
//     num_range = temp_ranges_thermo;
//     //plot all temperature ranges on one plot
//     for(int i=0; i<num_range; i++){

//       //specific polynomial
//       equation[i]<<R_UNIVERSAL/mol_mass<<"*("<<-thermo_data[i].Thermo_coef(0)<<"/(2.0*x*x) +"
// 		 <<-thermo_data[i].Thermo_coef(1)<<"/x +"
// 		 <<thermo_data[i].Thermo_coef(2)<<"*log(x)+"
// 		 <<thermo_data[i].Thermo_coef(3)<<"*x +"
// 		 <<thermo_data[i].Thermo_coef(4)<<"*x**2/2.0 +"
// 		 <<thermo_data[i].Thermo_coef(5)<<"*x**3/3.0 +"
// 		 <<thermo_data[i].Thermo_coef(6)<<"*x**4/4.0 +"
// 		 <<thermo_data[i].Thermo_coef(8)<<") "
// 		 <<" : 1/0 ";
      
//       //piecewise range using ternary opperator
//       segment[i]<<thermo_data[i].Low_range()<<"<=x && x<="<<thermo_data[i].High_range()<<" ? "; 
//     }
//     //overall range 
//     range<<" ["<<0<<":"<<thermo_data[num_range-1].High_range()<<"] "; 
//     units = "(J/(kg*K))"; 
//   }

//   //------------------ Viscosity -------------------------------------//
//   else if(type.compare("mu") == 0){
    
//     num_range = temp_ranges_V; 
//     //plot all temperature ranges on one plot
//     for(int i=0; i<num_range; i++){ 

//       //specific polynomial
//       equation[i]<<"exp("<<trans_viscosity[i].Trans_coef(0)<<"*log(x) +"
// 		 <<trans_viscosity[i].Trans_coef(1)<<"/x +" 
// 		 <<trans_viscosity[i].Trans_coef(2)<<"/(x*x) +"
// 		 <<trans_viscosity[i].Trans_coef(3)<<" )"  
// 		 <<" : 1/0 ";
      
//       //piecewise range using ternary opperator
//       segment[i]<<trans_viscosity[i].Low_range()<<"<=x && x<="<<trans_viscosity[i].High_range()<<" ? "; 
//     }
//     //overall range 
//     range<<" ["<<0<<":"<<trans_viscosity[num_range-1].High_range()<<"] "; 
//     units = "(kg/(m*s) or N*s/m^2)";
//   }

//   //--------------- Thermal Conductivity ----------------------------------//
//   else if(type.compare("k") == 0){
    
//     num_range = temp_ranges_C; 
//     //plot all temperature ranges on one plot
//     for(int i=0; i<num_range; i++){ 

//       //specific polynomial
//       equation[i]<<"exp("<<trans_thermconduct[i].Trans_coef(0)<<"*log(x) +"
// 		 <<trans_thermconduct[i].Trans_coef(1)<<"/x +" 
// 		 <<trans_thermconduct[i].Trans_coef(2)<<"/(x*x) +"
// 		 <<trans_thermconduct[i].Trans_coef(3)<<" )"  
// 		 <<" : 1/0 ";
      
//       //piecewise range using ternary opperator
//       segment[i]<<trans_thermconduct[i].Low_range()<<"<=x && x<="<<trans_thermconduct[i].High_range()<<" ? "; 
//     }
//     //overall range 
//     range<<" ["<<0<<":"<<trans_thermconduct[num_range-1].High_range()<<"] "; 
//     units = "(N/(s*K) or  W.(m*K))";
//   }

//   //---------------- ELSE -----------------------------------------------//
//   else{
//     cerr<<type<<" not a valid plot option (Cp,H,S,mu,k) \n";
//     exit(1);
//   }
  
//   //-----------------------------------------------------------------//
//   //--------- Create Script -----------------------------------------//
//   //-----------------------------------------------------------------//
//   //Create .gnu script file for 
//   gnu_out<<"set term post landscape color solid\n";
//   gnu_out<<"set out '"<<species<<"."<<type<<".ps'\n";
//   gnu_out<<"set title \"NASARP1311 "<<type<<" Plots for "<<species<<"\"\n";
//   gnu_out<<"set xlabel \"Temperature (K)\"\n";
//   gnu_out<<"set ylabel \""<<type<<" "<<units<<"\"\n";
//   gnu_out<<"set xzeroaxis\n";  
//   //using gnuplot ternary operator a?b:c to plot segments
//   //if a is true then b, else c
//   gnu_out<<"plot "<<range.str();
//   for(int i=0; i<num_range; i++){
//     gnu_out<<segment[i].str()<<equation[i].str();
//     if( i != num_range-1){
//       gnu_out<<",";
//     }
//   }
//   gnu_out<<"\n quit\n"; 
  
//   //Close file 
//   gnu_out.close(); 

//   //-----------------------------------------------------------------// 		
//   //--------- Run gnu plot scripts using system commands ------------//
//   //-----------------------------------------------------------------//
//   char buffer[80]; 
//   //run gnuplot scripts to create postscript plots
//   sprintf(buffer,"gnuplot %s",gnu_script);	   
//   system(buffer); 
   
//   //-------- Clean up files --------------------------------------/     
//   sprintf(buffer,"rm %s",gnu_script);  
//   system(buffer);
  
//   //clean up memory
//   delete[] equation;
//   delete[] segment;

// } //END Plot_data


/**************************************************************************
********************* LENNARD-JONES TRANPORT PROPS ************************
***************************************************************************/

// This is an addendum to the NASARP1311 species data class
// For additional information on the methods used to compute the transport 
// properties, see the Chemkin 4.0 Theory Manual (2004).


/************************************************************
 * Description: Viscosity computed Chapman-Enskog theory    * 
 *              See 'The properties of gases and liquids by *
 *              Reid et al.(1987).  This is Units [N s/m^2] *
 ************************************************************/
double NASARP1311data::Viscosity_LJ(double T) {
  
  /* Compute */
  double v = (5.0/16.0) * sqrt( PI * mol_mass  * BOLTZMANN * T / AVOGADRO ) 
    / ( PI * sigma*sigma * Omega22( T*BOLTZMANN/eps ) );
  
  return v;
}

// Polynomial fit to Viscosity_LJ
double NASARP1311data::Viscosity_LJ_Poly(double Temp){
  double logT = log(Temp);
  return exp( poly3(logT, trans_viscosity[0].trans_coef) );
}


double NASARP1311data::dViscositydT_LJ(double T) {
  
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
 ************************************************************/
double NASARP1311data::ThermalConduct_LJ(double T, double P) {
 
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
  double eta = Viscosity_LJ(T);
  
  
  /* Determine molar heat capacity relationships */
  
  Cv_trans = 3.0/2.0;
  
  // for a linear  molecule
  if (g==1) {
    Cv_rot = 1.0;
    Cv_vib = HeatCapacity_v(T)*mol_mass/R_UNIVERSAL - 5.0/2.0;
  }
  // for a non-linear  molecule
  else if (g==2) {
    Cv_rot = 3.0/2.0;
    Cv_vib = HeatCapacity_v(T)*mol_mass/R_UNIVERSAL - 3.0;  
  }

  /* Compute thermal conductivity */
  // for a monatomic gas
  if (g==0) {
    f_trans = 5.0/2.0;
    l = eta * R_UNIVERSAL * f_trans*Cv_trans / mol_mass;
  // for a non-monatonic
  } else {
    Dkk = BinaryDiff(*this,*this,T,P);  // self-diffusion coefficient
    v = eta;
    d = P*mol_mass / (R_UNIVERSAL*T);   // density
    a = d*Dkk/v;
    
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

// Polynomial fit to ThermalConduct_LJ
double NASARP1311data::ThermalConduct_LJ_Poly(double Temp, double Press){
  double logT = log(Temp);
  return exp( poly3(logT, trans_thermconduct[0].trans_coef) );
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
double Omega11(const double Ts) {

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
  if ( Ts<0.3 || Ts>100.0 ) {
    cerr << "NASARP1311dataclass::Omega11() - Warning: Reduced temperature T*= " 
	 << Ts << " outside applicable range (0.3 < T* < 100)." << endl;
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
double Omega22(const double Ts) {

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
  if ( Ts<0.3 || Ts>100.0 ) {
    cerr << "NASARP1311dataclass::Omega11() - Warning: Reduced temperature T*= " 
	 << Ts << " outside applicable range (0.3 < T* < 100)." << endl;
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
double dOmega22dT(const double T, const double eps) {
  
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
  if ( Ts<0.3 || Ts>100.0 ) {
    cerr << "Warning: Reduced temperature T*= " << Ts 
	 << " outside applicable range (0.3 < T* < 100)." << endl;
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
double F(const double Ts) {
  
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
double BinaryDiff( const NASARP1311data& s1, const NASARP1311data& s2, 
		   const double T, const double P) {
  
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
  const NASARP1311data* sp;    // polar species pointer
  const NASARP1311data* sn;    // non-polar species pointer
  
  
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
  double delta1, delta2;
  delta1 = 0.5*s1.mu*s1.mu/(s1.eps*s1.sigma*s1.sigma*s1.sigma);
  delta2 = 0.5*s2.mu*s2.mu/(s2.eps*s2.sigma*s2.sigma*s2.sigma);
  //cout << "\n delta1 = " << delta1 << "\t delta2 = " << delta2;


  // for both polar or both non-polar
  if ( ( fabs(delta1) > tol  &&  fabs(delta2) > tol)  ||  
       ( fabs(delta1) < tol  &&  fabs(delta2) < tol) ) {

    f_eps = f_sigma = 1.0;
      
    // for one polar and one non-polar
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

////////////////////////////////////////////////////////////////////////////////
///////////////////////FUNCTIONS////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/********************* String Converstions *************************/
// These are simple sting to ints and doubles conversion functions
// that work but are probably not the best way to do this.  
// string.h probably has a built in function to do just this.

int stoi(const string s)
{
  return atoi(s.c_str());
}

double stof(const string s)
{

  // strtod results in typeid errors using cxx "std:bad_typeid::what(void)"
  // if endp != endline character
  // also could use atof(s.c_str()), however no error checking
  char *endp;
  double Temp = strtod(s.c_str(),&endp); 

  //cout<<"\n "<<s<<" "<<s.c_str()<<" "<<Temp;

  if (s.c_str() != endp && *endp == '\0'){  
    //cout<<Temp<<endl;
    return Temp;
  } else {
    cerr<<"\n stof() error in reading NASARP1311dataclass::Getdata() "<<s; 
    exit(1);
  }
  
}




























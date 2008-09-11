/************************* Scalars.h **************************************

  This header file is for the Scalars class. Basically it is used to hold 
  the scalars that need to be solved in addition to the Favre-Filtered 
  Navier-Stokes equations.
               
*****************************************************************************/

#ifndef _SCALARS_INCLUDED   
#define _SCALARS_INCLUDED

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

/* Define some constants */
enum { NO_SCALARS, 
       LES_C, 
       LES_C_FSD, 
       LES_C_FSD_K, 
       LES_TF, 
       LES_TF_K };



class Set_scalar{

 private:
 public:
   string         *scalars;  //!< Scalars used in the turbulence model
   string    Scalar_system;  //!< Scalar system name
   static int  scalar_flag;  //!< Set flag of model scalars
   static int  num_scalars;  //!< Number of total scalars
  

   //! Constructor
   Set_scalar(){ scalars = NULL, Scalar_system = '\0'; }

   
   void scalar_set(string &);
   void Deallocate();


  //! Destructor
  ~Set_scalar(){ Deallocate(); };

};


/**************** Destructor *******************************************/
inline void Set_scalar::Deallocate(){
  //deallocate memory
  if(scalars != NULL) delete[] scalars; 
  scalars = NULL;
}



inline void Set_scalar::scalar_set(string &scal_sys) {

  Deallocate();

  Scalar_system = scal_sys;

  // Laminar case or DNS
  if( Scalar_system == "LAMINAR" || Scalar_system == "DNS"){
    scalar_flag = NO_SCALARS;
    num_scalars = 0;
    scalars = NULL;
  }

  // Progress variable model
  else if( Scalar_system == "LES_C"){
    scalar_flag = LES_C;
    num_scalars = 1;
    // set up scalar list
    scalars = new string[num_scalars];
    scalars[0] = "C";
  }

  // FSD & Progress variable model
  else if( Scalar_system == "LES_C_FSD"){
    scalar_flag = LES_C_FSD;
    num_scalars = 2;
    // set up scalar list
    scalars = new string[num_scalars];
    scalars[0] = "C";
    scalars[1] = "FSD/rho";
  }

  // FSD & Progress variable with SFS k model
  else if( Scalar_system == "LES_C_FSD_K"){
    scalar_flag = LES_C_FSD_K ;
    num_scalars = 3;
    // set up scalar list
    scalars = new string[num_scalars];
    scalars[0] = "C";
    scalars[1] = "FSD/rho";
    scalars[2] = "k";
  }

  // Thickened flame model
  else if( Scalar_system == "LES_TF"){
    scalar_flag = LES_TF;
    num_scalars = 0;
    // set up scalar list
    scalars = NULL;   
  }

  // Thickened flame with SFS k model
  else if( Scalar_system == "LES_TF_K"){
    scalar_flag = LES_TF_K;
    num_scalars = 1;
    // set up scalar list
    scalars = new string[num_scalars];
    scalars[0] = "K";
  }

  // else error  
  else{
    cerr<<"\n Set_scalar "<<Scalar_system<<" is not valid"<<endl;
    exit(1);
  }

} // end of scalar_set

#endif /* END of _SCALARS_INCLUDED */

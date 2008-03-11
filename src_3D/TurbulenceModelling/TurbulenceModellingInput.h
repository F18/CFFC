/* TurbulenceModellingInput.h:  Header file declaring turbulence model
                                input class. */

#ifndef _TURBULENCEMODEL_INPUT_INCLUDED
#define _TURBULENCEMODEL_INPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define TURBULENCEMODEL_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: Turbulence_Modelling_Input_Parameters
 *
 * @brief Input Parameters for the modelling of turbulence flows.
 *
 * This class defines and handles the input variables related to the
 * modelling and treatment of turbulent flows.
 *
 */
class Turbulence_Modelling_Input_Parameters{
 private:
 public:
    //@{ @name Indicator for near wall treatment for turbulent flows using FANS solver:
    //! Indicator for near wall treatment for turbulent flows (0, 1,2 , automatic, wall function, low_Reynolds number)
    int Near_Wall_Boundary_Treatment; 
    //! 0, 1,2 , automatic, wall function, low_Reynolds number.
    //@}
    
    //@{ @name LES related input parameters:
    char SFS_model[TURBULENCEMODEL_INPUT_PARAMETER_LENGTH];     //!< Sub-filter-scale model , default = Smagorinsky + Yoshizawa
    int i_SFS_model;                                            //!< Sub-filter-scale model , default = Smagorinsky + Yoshizawa
    char filter_type[TURBULENCEMODEL_INPUT_PARAMETER_LENGTH];   //!< Filter type : default = implicit filtering
    int i_filter_type;                                          //!< Filter type : default = implicit filtering
    double FGR;                                                 //!< Filter width to mesh size ratio, default : $f \bar{\Delta} = 2 \Delta x $f
    //@}
    
    //@{ @name Spectrum related input parameters:
    char spectrum[TURBULENCEMODEL_INPUT_PARAMETER_LENGTH];      //!< Spectrum function, default = Von Karman - Pao
    int i_spectrum;                                             //!< Spectrum function, default = Von Karman - Pao
    double LLR;                                                 //!< domainsize / (integral length scale L = k^(3/2)/eps)
    double TKE;                                                 //!< Turbulent Kinetic Energy
    //@}

  //@{ @name Constructors and desctructors:
  //! Constructor (assign default values)
  Turbulence_Modelling_Input_Parameters() {
      Near_Wall_Boundary_Treatment = 0; 
      // LES parameters
      strcpy(SFS_model, "Smagorinsky");
      i_SFS_model = SFS_MODEL_SMAGORINSKY;
      strcpy(filter_type, "Implicit");
      i_filter_type = FILTER_TYPE_IMPLICIT;
      FGR = TWO;
    
      // Spectrum parameters
      strcpy(spectrum,"Von_Karman_Pao");
      i_spectrum = SPECTRUM_VON_KARMAN_PAO;
      LLR = 6;
      TKE = 150.0;
  }

  //! Destructor
  ~Turbulence_Modelling_Input_Parameters(void){ }
  //@}

  //@{ @name Other Member functions:
  //! Broadcast input parameters to all processors
  void Broadcast(void);
  //! Parse next input line
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
  //! Check validity of specified input parameters
  int Check_Inputs(void);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
		               const Turbulence_Modelling_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Turbulence_Modelling_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  //@}

};

#endif // _TURBULENCEMODEL_INPUT_INCLUDED

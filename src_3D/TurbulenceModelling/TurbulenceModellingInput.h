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
    double smagorinsky_coefficient;                             //!< SFS smagorinsky coefficient
    int Filter_Initial_Condition;                               //!< Flag defines to filter the initial condition
    char filter_type[TURBULENCEMODEL_INPUT_PARAMETER_LENGTH];   //!< Filter type : default = implicit filtering
    int i_filter_type;                                          //!< Filter type : default = implicit filtering
    double FGR;                                                 //!< Filter width to mesh size ratio, default : $f \bar{\Delta} = 2 \Delta x $f
    double Filter_Width;                                        //!< Constant filter width
    int commutation_order;                                      //!< Commutation order of explicit filter
    int finite_differencing_order;                              //!< Finite differencing order in commutation error calculations
    int number_of_rings;                                        //!< Number of rings used in the explicit filter
    double Target_Filter_Sharpness;                             //!< Sharpness of target filter in case of Least squares constraint (vasilyev)
    bool Filter_Width_strict;                                   //!< This will strictly set the FGR and won't allow Least squares to approximate : default = false
    bool LS_constraints;                                        //!< This will turn on or off Least squares constraints for Vasilyev's filter : default = true
    int Derivative_constraints;                                 //!< Determines the number of derivative constraints for Vasilyev's filter : default = true
    bool Filter_Memory_Efficient;                               //!< Determines whether to store filter weights or not (not storing slows down tremendously)
    //@}
    
    //@{ @name Spectrum related input parameters:
    char spectrum[TURBULENCEMODEL_INPUT_PARAMETER_LENGTH];      //!< Spectrum function, default = Von Karman - Pao
    int i_spectrum;                                             //!< Spectrum function, default = Von Karman - Pao
    double LLR;                                                 //!< domainsize / (integral length scale L = k^(3/2)/eps)
    double TKE;                                                 //!< Turbulent Kinetic Energy
    int rescale_spectrum;                                       //!< Flag that decides if spectrum will be rescaled before computations
    //@}

    //@{ @name Reacting LES related input parameters:
    double           Laminar_Flame_Speed;                       //!< Propagation speed of laminar premixed flame
    double       Laminar_Flame_Thickness;                       //!< Thickness of laminar premixed flame
    double             Thickening_Factor;                       //!< Maximum thickening factor
    double        Fuel_Equivalence_Ratio;                       //!< Fuel equivalence ratio
    double    Unburnt_Fuel_Mass_Fraction;                       //!< Mass fraction of unburnt fuel (premixed flame)
    double             Reactants_Density;                       //!< Reactants Density
    double   Adiabatic_Flame_Temperature;                       //!< Adiabitc flame temperature
    //@}

  //@{ @name Constructors and desctructors:
  //! Constructor (assign default values)
  Turbulence_Modelling_Input_Parameters() {
      Near_Wall_Boundary_Treatment = 0; 
      // LES parameters
      strcpy(SFS_model, "Smagorinsky");
      i_SFS_model = SFS_MODEL_SMAGORINSKY;
      strcpy(filter_type, "Implicit");
      smagorinsky_coefficient = 0.18;
      i_filter_type = FILTER_TYPE_IMPLICIT;
      FGR = TWO;
      Filter_Width = ZERO;
      commutation_order = 2;
      finite_differencing_order = commutation_order+2;
      number_of_rings = 2;
      Target_Filter_Sharpness = -1;
      Filter_Initial_Condition = ON;
      Filter_Width_strict = OFF;
      LS_constraints = ON;
      Derivative_constraints = DEFAULT; // this lets an algorithm put the number
      Filter_Memory_Efficient = OFF;
    
      // Spectrum parameters
      strcpy(spectrum,"Pope");
      i_spectrum = SPECTRUM_POPE;
      LLR = 6;
      TKE = 150.0;
      rescale_spectrum = OFF;


      // Reacting LES parameters
      Laminar_Flame_Speed = 0.38;
      Laminar_Flame_Thickness = 0.44E-3;
      Thickening_Factor = ONE;
      Fuel_Equivalence_Ratio = ONE;
      Unburnt_Fuel_Mass_Fraction = 0.05518; 
      Reactants_Density = 1.13; 
      Adiabatic_Flame_Temperature = 2218.0;
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

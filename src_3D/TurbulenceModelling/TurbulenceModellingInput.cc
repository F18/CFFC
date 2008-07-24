/* TurbulenceModellingInput.cc: Definition of Turbulence_Modelling_Input_Parameters 
                                class member functions. */

/* Include the TurbulenceModellingInput header file. */

#ifndef _TURBULENCEMODEL_INPUT_INCLUDED
#include "TurbulenceModellingInput.h"
#endif // _TURBULENCEMODEL_INPUT_INCLUDED

/* Define member functions. */

/************************************************************************************
 * Turbulence_Modelling_Input_Parameters::Broadcast -- Broadcast to all processors. *
 ************************************************************************************/
void Turbulence_Modelling_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
    MPI::COMM_WORLD.Bcast(&(Near_Wall_Boundary_Treatment),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(SFS_model,
                          TURBULENCEMODEL_INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_SFS_model),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(smagorinsky_coefficient),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(filter_type,
                          TURBULENCEMODEL_INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_filter_type),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(FGR),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Width),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(commutation_order),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(number_of_rings),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Initial_Condition),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Target_Filter_Sharpness),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Width_strict),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(LS_constraints),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Derivative_constraints),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Filter_Memory_Efficient),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(spectrum,
                          TURBULENCEMODEL_INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_spectrum),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(LLR),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(TKE),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(rescale_spectrum),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Laminar_Flame_Speed),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Laminar_Flame_Thickness),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Thickening_Factor),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Fuel_Equivalence_Ratio),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Unburnt_Fuel_Mass_Fraction),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Reactants_Density),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Adiabatic_Flame_Temperature),
                          1,
                          MPI::DOUBLE, 0);

#endif

}

/********************************************************************************************
 * Turbulence_Modelling_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 ********************************************************************************************/
int Turbulence_Modelling_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                                              stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;

  if (strcmp(code, "Near_Wall_Boundary_Treatment") == 0) {
    i_command = 6001;
    value >> Near_Wall_Boundary_Treatment;
    if  (Near_Wall_Boundary_Treatment < 0 ||
	 Near_Wall_Boundary_Treatment > 2) i_command = INVALID_INPUT_VALUE;

    /*
     * LES parameters
     * --------------
     */
    
    /* ----- LES : Sub filter scale model ----- */
  } else if (strcmp(code, "SFS_model") == 0) {
    i_command = 120;
    value >> value_string;
    strcpy(SFS_model, value_string.c_str());
    if (strcmp(SFS_model, "Smagorinsky") == 0) {
        i_SFS_model = SFS_MODEL_SMAGORINSKY;
    } else if (strcmp(SFS_model, "k-equation") == 0) {
        i_SFS_model = SFS_MODEL_K_EQUATION;
    } else {
      i_command = INVALID_INPUT_VALUE;
    } 
      
  } else if (strcmp(code, "Smagorinsky_Coefficient") == 0) {
      i_command = 121;
      value >> smagorinsky_coefficient;
      if (smagorinsky_coefficient < ZERO)
          i_command = INVALID_INPUT_VALUE;
      
          
    /* ---- LES : filter type ---- */
  } else if (strcmp(code, "Filter_Type") == 0) {
      i_command = 130;
      value >> value_string;
      strcpy(filter_type, value_string.c_str());
      if (strcmp(filter_type, "Implicit") == 0) {
          i_filter_type = FILTER_TYPE_IMPLICIT;
      } else if (strcmp(filter_type, "Haselbacher") == 0) {
          i_filter_type = FILTER_TYPE_HASELBACHER;
      } else if (strcmp(filter_type, "Vasilyev") == 0) {
          i_filter_type = FILTER_TYPE_VASILYEV;
      } else if (strcmp(filter_type, "Tophat") == 0) {
          i_filter_type = FILTER_TYPE_TOPHAT;
      } else if (strcmp(filter_type, "Gaussian") == 0) {
          i_filter_type = FILTER_TYPE_GAUSSIAN;
      } else if (strcmp(filter_type, "Restart") == 0) {
          i_filter_type = FILTER_TYPE_RESTART;
      } else {
          i_command = INVALID_INPUT_VALUE;
      } /* endif */
      
  } else if (strcmp(code, "Filter_Memory_Efficient") == 0) {
      i_command = 132;
      value >> value_string;
      if (strcmp(value_string.c_str(), "ON") == 0) {
          Filter_Memory_Efficient = ON;
      } else if(strcmp(value_string.c_str(), "OFF") == 0) {
          Filter_Memory_Efficient = OFF; 
      } else {
          i_command = INVALID_INPUT_VALUE;
      }
      
  } else if (strcmp(code, "Filter_Grid_Ratio") == 0) {
    i_command = 131;
    value >> FGR;
    if (FGR < 1)
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Filter_Width") == 0) {
    i_command = 132;
    value >> Filter_Width;
    if ( Filter_Width < 0.0 )
      i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(code, "Target_Filter_Sharpness") == 0) {
      i_command = 132;
      value >> Target_Filter_Sharpness;
      
  } else if (strcmp(code, "Filter_Initial_Condition") == 0) {
      i_command = 132;
      value >> value_string;
      if (strcmp(value_string.c_str(), "ON") == 0) {
          Filter_Initial_Condition = ON;
      } else if(strcmp(value_string.c_str(), "OFF") == 0) {
          Filter_Initial_Condition = OFF; 
      } else {
          i_command = INVALID_INPUT_VALUE;
      }
      
  } else if (strcmp(code, "Filter_Width_strict") == 0) {
      i_command = 132;
      value >> value_string;
      if (strcmp(value_string.c_str(), "ON") == 0) {
          Filter_Width_strict = ON;
      } else if(strcmp(value_string.c_str(), "OFF") == 0) {
          Filter_Width_strict = OFF; 
      } else {
          i_command = INVALID_INPUT_VALUE;
      }
      
  } else if (strcmp(code, "LS_constraints") == 0) {
      i_command = 132;
      value >> value_string;
      if (strcmp(value_string.c_str(), "ON") == 0) {
          LS_constraints = ON;
      } else if(strcmp(value_string.c_str(), "OFF") == 0) {
          LS_constraints = OFF; 
      } else {
          i_command = INVALID_INPUT_VALUE;
      }
      
  } else if (strcmp(code, "Derivative_constraints") == 0) {
      i_command = 132;
      value >> Derivative_constraints;
      
  } else if (strcmp(code, "Commutation_Order") == 0) {
      i_command = 131;
      value >> commutation_order;
      if (commutation_order < 1)
          i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(code, "Number_of_Rings") == 0) {
      i_command = 131;
      value >> number_of_rings;
      if (number_of_rings < 1)
          i_command = INVALID_INPUT_VALUE;
        
    /* ---- Spectrum Parameters ---- */
  } else if (strcmp(code, "Spectrum_Model") == 0) {
    i_command = 140;
    value >> value_string;
    strcpy(spectrum, value_string.c_str());
    if (strcmp(spectrum, "Von_Karman_Pao") == 0) {
      i_spectrum = SPECTRUM_VON_KARMAN_PAO;
    } else if(strcmp(spectrum, "Haworth_Poinsot") == 0) {
      i_spectrum = SPECTRUM_HAWORTH_POINSOT;
    } else if(strcmp(spectrum, "Pope") == 0) {
      i_spectrum = SPECTRUM_POPE;
    } else if(strcmp(spectrum, "Laval_Nazarenko") == 0) {
      i_spectrum = SPECTRUM_LAVAL_NAZARENKO; 
    } else if(strcmp(spectrum, "Uniform") == 0) {
      i_spectrum = SPECTRUM_UNIFORM;
    } else {
      i_command = INVALID_INPUT_VALUE;
    } /* endif */
        
  } else if (strcmp(code, "Domain_Integral_Lengthscale_Ratio") == 0) {
    i_command = 141;
    value >> LLR;
    if (LLR < 1)
      i_command = INVALID_INPUT_VALUE;
        
  } else if (strcmp(code, "Turbulent_Kinetic_Energy") == 0) {
    i_command = 142;
    value >> TKE;
    if (TKE <= ZERO)
      i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(code, "Rescale_Spectrum_Before_Execution") == 0) {
      i_command = 143;
      value >> value_string;
      if      (value_string == "ON")   rescale_spectrum = ON;
      else if (value_string == "OFF")  rescale_spectrum = OFF;
      else 
          i_command = INVALID_INPUT_VALUE;

    /* ----- Reacting LES parameters ----- */
  } else if (strcmp(code, "Laminar_Flame_Speed") == 0) {
    i_command = 150;
    value >> Laminar_Flame_Speed;
    if ( Laminar_Flame_Speed <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Laminar_Flame_Thickness") == 0) {
    i_command = 151;
    value >> Laminar_Flame_Thickness;
    if ( Laminar_Flame_Thickness <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Thickening_Factor") == 0) {
    i_command = 152;
    value >> Thickening_Factor;
    if ( Thickening_Factor < ONE )
      i_command = INVALID_INPUT_VALUE;
    
  } else if (strcmp(code, "Fuel_Equivalence_Ratio") == 0) {
    i_command = 153;
    value >> Fuel_Equivalence_Ratio;
    if ( Fuel_Equivalence_Ratio <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Unburnt_Fuel_Mass_Fraction") == 0) {
    i_command = 154;
    value >> Unburnt_Fuel_Mass_Fraction;
    if ( Unburnt_Fuel_Mass_Fraction <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Reactants_Density") == 0) {
    i_command = 155;
    value >> Reactants_Density;
    if ( Reactants_Density <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(code, "Adiabatic_Flame_Temperature") == 0) {
    i_command = 156;
    value >> Adiabatic_Flame_Temperature;
    if ( Adiabatic_Flame_Temperature <= 0.0 )
      i_command = INVALID_INPUT_VALUE;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/******************************************************************************
 * Turbulence_Modelling_Input_Parameters::Check_Inputs -- Check input values. *
 ******************************************************************************/
int Turbulence_Modelling_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * Turbulence_Modelling_Input_Parameters -- Input-output operators.        *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Turbulence_Modelling_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      Turbulence_Modelling_Input_Parameters &IP) {

  return in_file;

}

void Turbulence_Modelling_Input_Parameters::Output(ostream &out_file) const {

  if (Near_Wall_Boundary_Treatment == 0) {
     out_file << "\n  -> Near Wall Turbulent BC Treatment: Automatic";
  } else {
     if (Near_Wall_Boundary_Treatment == 2) {
        out_file << "\n  -> Near Wall Turbulent BC Treatment: Direct Integration";
     } else {
        out_file << "\n  -> Near Wall Turbulent BC Treatment: Wall Functions";
     } /* endif */
  } /* endif */

    out_file << "\n  LES parameters:";
    out_file << "\n    -> Filter type: " << filter_type;
    if (i_filter_type == FILTER_TYPE_IMPLICIT){
        out_file << "\n       -> Filter Grid Ratio: " << FGR;
    }
    out_file << "\n    -> Sub Filter Scale model: " << SFS_model;
    if (i_SFS_model == SFS_MODEL_SMAGORINSKY){
        out_file << "\n       -> Smagorinsky Coefficient: " << smagorinsky_coefficient;
    }
    
    out_file << "\n  Spectrum parameters:";
    out_file << "\n    -> Spectrum model: " << spectrum;
    out_file << "\n    -> Domain Integral length scale \"L\" ratio: " << LLR;
    out_file << "\n    -> Turbulent Kinetic Energy: " << TKE; 

}

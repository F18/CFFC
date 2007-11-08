/*!\file DiffusionFields.h
  \brief Header file defining 2D diffusion coefficient fields. */

#ifndef _DIFFUSION_FIELDS_INCLUDED
#define _DIFFUSION_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../CFD/CFD.h"

/*!
 * \class DiffusionFields
 *
 * @brief Collection of diffusion coefficient fields for 2D problems
 *
 * \nosubgrouping
 */

class DiffusionFields{
public:
  //! @name Defined public types
  //@{
  //! Non-linear diffusion field type
  typedef double (*NonlinearDiffusionFieldType)(const double &, const double &, const double &);

  //! Reference to non-linear diffusion field type
  typedef double (*&NonlinearDiffusionFieldTypeRef)(const double &, const double &, const double &);
  //@}

  //! @name Defined non-linear diffusion fields
  //@{
  static double Zero_Diffusion(const double &x, const double &y, const double &u);
  static double Constant_Diffusion(const double &x, const double &y, const double &u);
  static double Linear_Diffusion(const double &x, const double &y, const double &u);
  //@}

  //! @name Functions that set class parameters
  //@{
  static void Set_ConstantDiffusionField(const double & DiffusionCoeff);
  //@}

  //! @name Functions for input-output and broadcast
  //@{
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);
  static void Print_Info(std::ostream & out_file);
  static void Broadcast(void);
  //@}

  //! @name Functions for setting external pointers
  //@{
  static void Connect_Pointer_To_Diffusion_Field(NonlinearDiffusionFieldTypeRef DiffusionField);
  //@}

protected:
  DiffusionFields(void);   //!< Private default constructor
  DiffusionFields(const DiffusionFields&); //!< Private copy constructor
  DiffusionFields& operator=(const DiffusionFields&); //!< Private assignment operator

  //! @name Parameters used to determine a particular type of diffusion field
  //@{
  static short i_Diffusion_Field_Type;     //!< type indicator for the diffusion coefficient field
  static double DiffusionCoefficient;      //!< reference diffusion coefficient
  //@}

  //! @name Set field types
  //@{
  static int Set_Diffusion_Field_Type(const char * FieldType);
  //@}

};

/*! Defines a ZERO diffusion field */
inline double DiffusionFields::Zero_Diffusion(const double &x, const double &y, const double &u){
  return 0.0;
}

/*! Defines a constant diffusion field */
inline double DiffusionFields::Constant_Diffusion(const double &x, const double &y, const double &u){
  return DiffusionCoefficient;
}

/*!
 * Defines a linear diffusion field (varies linearly in both x-direction and y-direction)
 *
 * \todo Implement this field as: k(x,y) = k_x*(x-xr) + k_y*(y-yr) + kr
 */
inline double DiffusionFields::Linear_Diffusion(const double &x, const double &y, const double &u){
  return 0.0;
}

//! Parse the input control parameters for 
//  settings related to DiffusionFields class
template<class Input_Parameters_Type> inline
void DiffusionFields::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Diffusion_Field") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    i_command = Set_Diffusion_Field_Type(IP.Next_Control_Parameter);
  } else if (strcmp(IP.Next_Control_Parameter, "Diffusion_Coefficient") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> DiffusionCoefficient;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif

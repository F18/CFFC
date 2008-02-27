/*!\file AdvectDiffuse2DInflowField.h
  \brief Header file defining the singleton advection-diffusion exact solution class. */

#ifndef _ADVECTDIFFUSE2D_INFLOWFIELD_SINGLETON_INCLUDED
#define _ADVECTDIFFUSE2D_INFLOWFIELD_SINGLETON_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "InflowFields.h"

/*! 
 * \class AdvectDiffuse2D_InflowField
 * 
 * \brief Collection of inflow fields for Advection-Diffusion 2D problems
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class AdvectDiffuse2D_InflowField{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static AdvectDiffuse2D_InflowField& getInstance(void);
  void SetInflowField(const short &InflowIndex);
  static void DestroyInflowFieldObject(void);
  //@}

  //! @name Evaluation of exact solution
  //@{ 
  double operator()(const double &x, const double &y) const;
  double Solution(const double &x, const double &y) const;
  //@}

  //! Update the internal variables of the exact solution
  void Set_InflowField_Parameters(void) const {
    if ( IsInflowFieldSet() ){ return Inflow->Set_InflowField_Parameters(); } 
  }

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

  //! Indicate whether the inflow field has been set or not
  bool IsInflowFieldSet(void) const { return (Inflow != NULL)? true: false; }

protected:
  AdvectDiffuse2D_InflowField();		//!< Private constructor
  AdvectDiffuse2D_InflowField(const AdvectDiffuse2D_InflowField&); //!< Private copy constructor
  AdvectDiffuse2D_InflowField& operator=(const AdvectDiffuse2D_InflowField&); //!< Private assignment operator
  
private:
  static InflowFieldBasicType *Inflow;    //!< pointer to the inflow field
  static short i_Inflow_Field_Type;       //!< type indicator for the inflow field
};

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_InflowField::operator()(const double &x, const double &y) const{
  return Inflow->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_InflowField::Solution(const double &x, const double &y) const{
  return Inflow->EvaluateSolutionAt(x,y);
}

#endif

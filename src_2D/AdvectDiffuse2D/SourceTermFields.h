/*!\file SourceTermFields.h
  \brief Header file defining 2D source term fields. */

#ifndef _SOURCETERM_FIELDS_INCLUDED
#define _SOURCETERM_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "SourceFields.h"

/*! 
 * \class SourceTermFields
 * 
 * \brief Collection of source term fields for Advection-Diffusion 2D problems
 *
 * The source term fields are considered to have one of the following forms:\n
 *
 * - Pointwise form: \f$ \phi= \phi(x,y,u(x,y)) \f$ , x,y Cartesian coordinates, u(x,y) solution
 * - Integral form as function of the average solution over the computed domain, \f$\bar{u}\f$, : \f$ \phi= \phi(\bar{u}) \f$
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class SourceTermFields{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static SourceTermFields& getInstance(void);
  void SetSourceField(const short &FieldIndex);
  static void DestroySource(void);
  //@}

  //! @name Evaluation of field solution
  //@{ 
  double operator()(const double &x, const double &y, const double &u);
  double operator()(const double &u_avg);
  //@}

  //! Specifies if the current field require or not integration
  bool FieldRequireIntegration(void){ return SourcePtr->FieldRequireIntegration();}

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

protected:
  SourceTermFields();		//!< Private constructor
  SourceTermFields(const SourceTermFields&); //!< Private copy constructor
  SourceTermFields& operator=(const SourceTermFields&); //!< Private assignment operator
  
private:
  static SourceFieldBasicType *SourcePtr; //!< pointer to the source field
  static short i_Source_Field_Type;       //!< type indicator for the source field
};

/*! 
 * Evaluate the source term field for the pointwise formulation
 * \param u the solution value at the x and y Cartesian coordinates
 */
inline double SourceTermFields::operator()(const double &x, const double &y, const double &u){
  return SourcePtr->FieldSoln(x,y,u);
}

/*! 
 * Evaluate the source term field for the integral formulation
 * \param u_avg the average solution
 */
inline double SourceTermFields::operator()(const double &u_avg){
  return SourcePtr->FieldSoln(u_avg);
}

#endif

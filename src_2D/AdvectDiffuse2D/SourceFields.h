/*!\file SourceFields.h
  \brief Header file defining 2D source fields. */

#ifndef _SOURCE_FIELDS_INCLUDED
#define _SOURCE_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Math.h"
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

// Declare the input parameters class
class AdvectDiffuse2D_Input_Parameters;

/*! 
 * \class SourceFieldBasicType
 *
 * \brief Basic data type for any source field.
 *
 * This is an abstract data type (ADT).
 */
class SourceFieldBasicType{
public:

  //! Default ctor
  SourceFieldBasicType(void);

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~SourceFieldBasicType(void) = 0;

  /*! Calculate the non-linear field solution based on the location and the solution value */
  virtual double FieldSoln(const double &x, const double &y, const double & Soln) = 0;

  /*! Calculate the non-linear field solution based on the average solution value */
  virtual double FieldSoln(const double &SolnAvg) = 0;

  /*! Calculate the stability limit for the non-linear source field */
  virtual double StabilityLimit(const double &x, const double &y, const double & Soln) = 0;

  /*! Calculate the stability limit for the non-linear source field */
  virtual double StabilityLimit(const double &SolnAvg) = 0;

  /*! This function returns true is the solution needs to be integrated over the cell domain and false if it does not. */
  virtual bool FieldRequireIntegration(void) = 0;

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command) = 0;
  virtual void Print_Info(std::ostream & out_file) = 0;
  virtual void Broadcast(void) = 0;
  //@}
  //@}

  string & whatField(void){ return FieldName; }   //!< Get the name of the source field type

protected:
  string FieldName;		//!< storage for the name of the source field type

};

/*! 
 * \class SourceFieldThatRequireIntegration
 *
 * \brief Basic data type for pointwise source fields.
 *
 * This type of field requires numerical integration of the field solution
 * in order to obtain the average source solution over a domain.
 */
class SourceFieldThatRequireIntegration: public SourceFieldBasicType{
public:

  //! This field require integration (information used at higher level)
  bool FieldRequireIntegration(void){ return true; }

  //! Declare a pure virtual destructor
  virtual ~SourceFieldThatRequireIntegration(void) = 0;

  //! \throw runtime_error if this function gets called for this class
  double FieldSoln(const double &SolnAvg){
    throw runtime_error("SourceFieldThatRequireIntegration::FieldSoln() ERROR! Wrong number of parameters!");
    return 0.0;
  }

  //! \throw runtime_error if this function gets called for this class
  double StabilityLimit(const double &SolnAvg){
    throw runtime_error("SourceFieldThatRequireIntegration::StabilityLimit() ERROR! Wrong number of parameters!");
    return 0.0;
  }
};

/*! 
 * \class SourceFieldThatDoesNotRequireIntegration
 *
 * \brief Basic data type for integral source fields.
 *
 * This type of field DON'T require numerical integration of the field solution
 * in order to obtain the average source solution over a domain.
 * The integral of the field solution is based on the average solution and it also 
 * includes the area dependency.
 */
class SourceFieldThatDoesNotRequireIntegration: public SourceFieldBasicType{
public:

  //! This field DOES NOT require integration (information used at higher level)
  bool FieldRequireIntegration(void){ return false; }

  //! Declare a pure virtual destructor
  virtual ~SourceFieldThatDoesNotRequireIntegration(void) = 0;

  //! \throw runtime_error if this function gets called for this class
  double FieldSoln(const double &x, const double &y, const double & Soln){
    throw runtime_error("SourceFieldThatDoesNotRequireIntegration::FieldSoln() ERROR! Wrong number of parameters!");
    return 0.0;
  }

  //! \throw runtime_error if this function gets called for this class
  double StabilityLimit(const double &x, const double &y, const double & Soln){
    throw runtime_error("SourceFieldThatRequireIntegration::StabilityLimit() ERROR! Wrong number of parameters!");
    return 0.0;
  }
};


/*! 
 * \class ZERO_SourceField
 * 
 * \brief Implements a ZERO value source field (no source).
 *
 * The field has the expression: \f$ \phi= \phi(\bar{u}) = 0 \f$
 */
class ZERO_SourceField: public SourceFieldThatDoesNotRequireIntegration{
public:

  //! Basic Constructor
  ZERO_SourceField(void){ 
    FieldName = "No source";	// Name the field
  }

  //! Return ZERO solution (no source)
  double FieldSoln(const double &SolnAvg){return 0.0;}

  //! Return the maximum allowed time step which ensures a stable computation for this field (i.e. stability limit)
  double StabilityLimit(const double &SolnAvg){ return MILLION;}

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command){};

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file){};

  //! Broadcast relevant parameters
  void Broadcast(void){};
};

/*! 
 * \class Linear_SourceField
 * 
 * \brief Implements a linear source field with a relaxation coefficient
 *         
 * The field has the expression: \f$ \phi= \phi(\bar{u}) = \bar{u}/\tau \f$
 */
class Linear_SourceField: public SourceFieldThatDoesNotRequireIntegration{
public:

  //! Basic Constructor
  Linear_SourceField(void): tau(1.0){
    FieldName = "Linear variation [u/tau]";	// Name the field
  };

  //! Return linear solution 
  double FieldSoln(const double &SolnAvg);

  //! Return the maximum allowed time step which ensures a stable computation for this field (i.e. stability limit)
  double StabilityLimit(const double &SolnAvg){ return HALF*tau;}

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double tau;			//!< relaxation coefficient
};

/*!
 * Implement the source field expression
 */
inline double Linear_SourceField::FieldSoln(const double &SolnAvg){
  return SolnAvg/tau;
}


/*! 
 * \class Exponential_SourceField
 * 
 * \brief Implements a non-linear exponential source field.
 * 
 * The field has the expression: \f$ \phi= \phi(x,y,u(x,y)) = a\, e^{\beta \, u} \f$
 */
class Exponential_SourceField: public SourceFieldThatRequireIntegration{
public:
  
  Exponential_SourceField(void): a(1.0), beta(1.0){
    FieldName = "Exponential variation [A*e^(B*u)]";	// Name the field
  };

  //! Return exponential non-linear solution
  double FieldSoln(const double &x, const double &y, const double & Soln);

  //! Return the maximum allowed time step which ensures a stable computation for this field (i.e. stability limit)
  double StabilityLimit(const double &x, const double &y, const double & Soln);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double a, beta;		//!< coefficients of the field
};

/*!
 * Implement the source field expression
 */
inline double Exponential_SourceField::FieldSoln(const double &x, const double &y, const double & Soln){
  return a*exp(beta*Soln);
}

/*!
 * Return the maximum allowed time step which ensures a stable computation for this field (i.e. stability limit)
 *
 * \todo Work out the stability limit
 */
inline double Exponential_SourceField::StabilityLimit(const double &x, const double &y, const double & Soln){
  return MILLION;
}

#endif

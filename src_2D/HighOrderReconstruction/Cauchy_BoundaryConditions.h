/*!\file Cauchy_BoundaryConditions.h
  \brief Header file defining Cauchy_BCs template class for different solution types. */

/*************************************************************************************************************
  Cauchy type boundary conditions specify a weighted average of Dirichlet and Neumann kinds.
  Therefore this template allocates memory for the values of the function on a surface (T=f(r,t)) as well
  as for the normal derivatives of the function on the same surface.
  The boundary conditions can be specified at a variable number of locations.
*************************************************************************************************************/

#ifndef _CAUCHY_BOUNDARYCONDITIONS_INCLUDED
#define _CAUCHY_BOUNDARYCONDITIONS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Utilities/Utilities.h"


template< class SOLN_STATE >
class Cauchy_BCs;

/************************************************
*     Friend Functions :   Cauchy_BCs          *
************************************************/
template<class SOLN_STATE >
std::ostream& operator<< (std::ostream& os, const Cauchy_BCs<SOLN_STATE> & Obj);

template<class SOLN_STATE >
std::istream& operator>> (std::istream& is, Cauchy_BCs<SOLN_STATE> & Obj);

template< class SOLN_STATE >
class Cauchy_BCs{

public:

  //! @name Defined public types:
  //@{
  typedef SOLN_STATE Soln_State;
  //@}

  //! @name Constructors:
  //@{
  Cauchy_BCs(void);
  Cauchy_BCs(const Cauchy_BCs<SOLN_STATE> & rhs);
  //@}

  //! @name Destructors:
  //@{
  ~Cauchy_BCs(void){ deallocate(); }
  void deallocate(void);
  //@}

  Cauchy_BCs<SOLN_STATE> & operator=(const Cauchy_BCs<SOLN_STATE>& rhs); //!< Assignment operator

  //! @name Initialize container functions.
  //@{
  void allocate(const int & NumOfObjects);
  void InitializeCauchyBCs(const int & NumOfObjects,
			   const int & BCtype);
  //@}

  //! @name Set number of constraints
  //@{
  //! @brief Set number of constraints corresponding to the type of boundary condition
  void SetCharacteristicConstraintsBasedOnBCtype(const int & BCtype);
  //@}

  //! @name Field access:
  //@{
  const int & NumOfPoints(void) const {return NumLoc;}

  const SOLN_STATE * DirichletBC(void){return Dirichlet; }
  SOLN_STATE & DirichletBC(const int &position){return Dirichlet[position-1];}
  const SOLN_STATE & DirichletBC(const int &position) const {return Dirichlet[position-1];}

  const SOLN_STATE * NeumannBC(void) {return Neumann; }
  SOLN_STATE & NeumannBC(const int &position) {return Neumann[position-1];}
  const SOLN_STATE & NeumannBC(const int &position) const {return Neumann[position-1];}

  const SOLN_STATE * a(void){return a_coeff; }
  SOLN_STATE & a(const int &position){return a_coeff[position-1];}
  const SOLN_STATE & a(const int &position) const {return a_coeff[position-1];}

  const SOLN_STATE * b(void) {return b_coeff; }
  SOLN_STATE & b(const int &position) {return b_coeff[position-1];}
  const SOLN_STATE & b(const int &position) const {return b_coeff[position-1];}
  
  const ShortIndexType & IndividualConstraintsVector(void) const { return IndividualConstraints; }
  short int & NumberOfIndividualConstraints(const int &parameter){ return IndividualConstraints[parameter-1]; }
  const short int & NumberOfIndividualConstraints(const int &parameter) const { return IndividualConstraints[parameter-1]; }

  const ShortIndexType & RelationalConstraintsVector(void) const { return RelationalConstraints; }
  short int & NumberOfRelationalConstraints(const int &parameter){ return RelationalConstraints[parameter-1]; }
  const short int & NumberOfRelationalConstraints(const int &parameter) const { return RelationalConstraints[parameter-1]; }

  bool IsThereAnyRelationalConstraintRequired(void) const { return RelationalConstraintsFlag; }
  bool IsThereAnyIndividualConstraintRequired(void) const { return IndividualConstraintsFlag; }

  void setIndividualConstraints(void) { IndividualConstraintsFlag = true; }
  void setNOIndividualConstraints(void) { IndividualConstraintsFlag = false; }
  void setRelationalConstraints(void) { RelationalConstraintsFlag = true; }
  void setNORelationalConstraints(void) { RelationalConstraintsFlag = false; }
  

  //! @name Friend functions:
  //@{
  friend std::ostream& operator<< <SOLN_STATE> (std::ostream& os, const Cauchy_BCs<SOLN_STATE>& rhs);
  friend std::istream& operator>> <SOLN_STATE> (std::istream& is, Cauchy_BCs<SOLN_STATE>& rhs);
  //@}

private:
  SOLN_STATE *Dirichlet;	//!< Dirichlet boundary conditions
  SOLN_STATE *Neumann;	        //!< Neumann boundary conditions
  SOLN_STATE *a_coeff;		//!< Coefficient of Dirichlet BC in the mixed BC
  SOLN_STATE *b_coeff;		//!< Coefficient of Neumann BC in the mixed BC
  int NumLoc;			//!< Number of locations where the BCs are specified

  //! @name Vectors of constraints for each solution parameter at each particular Location:
  //@{  
  ShortIndexType IndividualConstraints; /*!< Number of individual constraints imposed by the type of BCs.
					 *   Get the total number of individual constraints by multiplying with NumLoc. */
  ShortIndexType RelationalConstraints; /*!< Number of relational constraints imposed by the type of BCs.
					 *   Get the total number of individual constraints by multiplying with NumLoc. */
  //@}

  bool IndividualConstraintsFlag; //!< Flag for existence of individual constraints
  bool RelationalConstraintsFlag; //!< Flag for existence of relational constraints

  bool _allocated_container;	//!< allocation flag
};

/*******************************************************
 * CLASS Template:   Cauchy_BCs                        *
 * Implementation of the Member Functions              *
 ******************************************************/

//! Default Constructor
template< class SOLN_STATE > inline
Cauchy_BCs<SOLN_STATE>::Cauchy_BCs(void):
  Dirichlet(NULL), Neumann(NULL),
  a_coeff(NULL), b_coeff(NULL),
  NumLoc(0),
  IndividualConstraintsFlag(false),
  RelationalConstraintsFlag(false),
  _allocated_container(false)
{
  //
}

//! Copy constructor 
template< class SOLN_STATE > inline
Cauchy_BCs<SOLN_STATE>::Cauchy_BCs(const Cauchy_BCs<SOLN_STATE> & rhs):
  Dirichlet(NULL), Neumann(NULL),
  a_coeff(NULL), b_coeff(NULL),
  NumLoc(0),
  IndividualConstraintsFlag(false),
  RelationalConstraintsFlag(false)
{

  // check if the rhs has memory allocated
  if (rhs._allocated_container){

    /* allocate memory for the new object */
    allocate(rhs.NumOfPoints());

    /* copy the values from the RHS */
    for (int i=0; i<rhs.NumOfPoints(); ++i){
      Dirichlet[i] = rhs.Dirichlet[i];
      Neumann[i]   = rhs.Neumann[i];
      a_coeff[i]   = rhs.a_coeff[i];
      b_coeff[i]   = rhs.b_coeff[i];      
    }

    /* copy constraints values from the RHS */
    IndividualConstraints = rhs.IndividualConstraints;
    RelationalConstraints = rhs.RelationalConstraints;

    IndividualConstraintsFlag = rhs.IndividualConstraintsFlag;
    RelationalConstraintsFlag = rhs.RelationalConstraintsFlag;
  }
}

/*!
 * Allocate memory for the object.
 */
template< class SOLN_STATE > inline
void Cauchy_BCs<SOLN_STATE>::allocate(const int &NumOfObjects){

  if (NumLoc == NumOfObjects) {
    return; /* enough memory already allocated */
  } 

  /* If there is not enough memory, deallocate the current one and allocate again */
  if (_allocated_container){
    deallocate();
  }

  // allocate new memory
  if (NumOfObjects > 0){	
    NumLoc = NumOfObjects;
    Dirichlet = new SOLN_STATE[NumLoc];
    Neumann   = new SOLN_STATE[NumLoc];
    a_coeff   = new SOLN_STATE[NumLoc];
    b_coeff   = new SOLN_STATE[NumLoc];

    // allocate memory for the constraints vectors
    IndividualConstraints.reserve(SOLN_STATE::NumVar());
    RelationalConstraints.reserve(SOLN_STATE::NumVar());

    // Initialize constraints vectors
    for (int parameter = 1; parameter <= SOLN_STATE::NumVar(); ++parameter){
      IndividualConstraints.push_back(0); // no individual constraints for this parameter
      RelationalConstraints.push_back(0); // no relational constraints for this parameter
    }
  }
  
  // Confirm the allocation
  _allocated_container = true;
}

/*!
 * Initialize container with a particular boundary condition
 * and a specified number of objects.
 */
template< class SOLN_STATE > inline
void Cauchy_BCs<SOLN_STATE>::InitializeCauchyBCs(const int &NumOfObjects,
						 const int &BCtype){

  // allocate memory
  allocate(NumOfObjects);

  if (NumOfObjects > 0){
    // set number of constraints per GQP for each solution parameter
    SetCharacteristicConstraintsBasedOnBCtype(BCtype);
  }
}


/*!
 * Deallocate memory
 */
template< class SOLN_STATE > inline
void Cauchy_BCs<SOLN_STATE>::deallocate(void){
  delete [] Dirichlet; Dirichlet = NULL;
  delete [] Neumann; Neumann = NULL;
  delete [] a_coeff; a_coeff = NULL;
  delete [] b_coeff; b_coeff = NULL;
  NumLoc = 0;

  IndividualConstraints.clear();
  RelationalConstraints.clear();

  IndividualConstraintsFlag = false;
  RelationalConstraintsFlag = false;

  // Confirm the de-allocation
  _allocated_container = false;
}

/*!
 * Assignment operator
 */
template< class SOLN_STATE > inline
Cauchy_BCs<SOLN_STATE> & Cauchy_BCs<SOLN_STATE>::operator=(const Cauchy_BCs<SOLN_STATE>& rhs){

  // !!! If the LHS container already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & rhs) return *this;

  // allocate memory
  allocate(rhs.NumOfPoints());

  /* copy the values from the RHS */
  for (int i=0; i<rhs.NumOfPoints(); ++i){
    Dirichlet[i] = rhs.Dirichlet[i];
    Neumann[i]   = rhs.Neumann[i];
    a_coeff[i]   = rhs.a_coeff[i];
    b_coeff[i]   = rhs.b_coeff[i];      
  }

  /* copy constraints values from the RHS */
  IndividualConstraints = rhs.IndividualConstraints;
  RelationalConstraints = rhs.RelationalConstraints;
  
  IndividualConstraintsFlag = rhs.IndividualConstraintsFlag;
  RelationalConstraintsFlag = rhs.RelationalConstraintsFlag;

  return *this;
}

/*!
 * Set number of constraints (i.e. individual and relational) 
 * required by the boundary condition type to be satisfied  
 * on each solution parameter at each flux calculation point (GQP).
 *
 * \param BCtype specifies the boundary condition type
 *
 * \note Specialize this routine to set the proper number
 *       of constraints for each solution state!
 */
template< class SOLN_STATE > inline
void Cauchy_BCs<SOLN_STATE>::SetCharacteristicConstraintsBasedOnBCtype(const int & BCtype){

  for (int parameter = 1; parameter <= IndividualConstraints.size(); ++parameter) {

    // Impose NO relational constraints
    RelationalConstraints[parameter-1] = 0;

    // Impose No individual constraints
    IndividualConstraints[parameter-1] = 0;
  }

  IndividualConstraintsFlag = false;
  RelationalConstraintsFlag = false;
}

/* Friend functions */

//! operator<<
template< class SOLN_STATE > inline
std::ostream& operator<< (std::ostream& os, const Cauchy_BCs<SOLN_STATE>& rhs){

  os << rhs.NumOfPoints() << std::endl;
  os.precision(15);
  for (int i=1; i<=rhs.NumOfPoints(); ++i){
    os << rhs.DirichletBC(i) << "\n"
       << rhs.NeumannBC(i)   << "\n"
       << rhs.a(i) << "\n"
       << rhs.b(i) << "\n";
  }

  return os;
}

//! operator>>
template< class SOLN_STATE > inline
std::istream& operator>> (std::istream& is, Cauchy_BCs<SOLN_STATE>& rhs){

  int N(0);

  is.setf(ios::skipws);
  /* Read the number of points */
  is >> N;

  /* Allocate memory if the current memory has a different dimension */
  rhs.allocate(N);

  /* Read the data */
  for (int i=1; i<=N; ++i){
    is >> rhs.DirichletBC(i)
       >> rhs.NeumannBC(i)
       >> rhs.a(i)
       >> rhs.b(i);
  }

  return is;
}


/*****************************************************************
 *            SPECIALIZATIONS FOR CLASS DOUBLE                   *
 ****************************************************************/

/*!
 * Allocate memory for the object.
 */
template< > inline
void Cauchy_BCs<double>::allocate(const int &NumOfObjects){

  if (NumLoc == NumOfObjects) {
    return; /* enough memory already allocated */
  } 

  /* If there is not enough memory, deallocate the current one and allocate again */
  if (_allocated_container){
    deallocate();
  }

  // allocate new memory
  if (NumOfObjects > 0){
    NumLoc = NumOfObjects;
    Dirichlet = new double[NumLoc];
    Neumann   = new double[NumLoc];
    a_coeff   = new double[NumLoc];
    b_coeff   = new double[NumLoc];

    // allocate memory for the constraints vectors
    IndividualConstraints.reserve(1);
    RelationalConstraints.reserve(1);

    // Initialize constraints vectors
    IndividualConstraints.push_back(0);
    RelationalConstraints.push_back(0);
  }
  
  // Confirm the allocation
  _allocated_container = true;
}

/*!
 * Access individual constraints
 */
template< > inline
short int & Cauchy_BCs<double>::NumberOfIndividualConstraints(const int &parameter){
  return IndividualConstraints[0];
}

/*!
 * Access individual constraints with 'const'
 */
template< > inline
const short int & Cauchy_BCs<double>::NumberOfIndividualConstraints(const int &parameter) const{
  return IndividualConstraints[0];
}

/*!
 * Access relational constraints
 */
template< > inline
short int & Cauchy_BCs<double>::NumberOfRelationalConstraints(const int &parameter){
  return RelationalConstraints[0];
}

/*!
 * Access relational constraints with 'const'
 */
template< > inline
const short int & Cauchy_BCs<double>::NumberOfRelationalConstraints(const int &parameter) const{
  return RelationalConstraints[0];
}

#endif

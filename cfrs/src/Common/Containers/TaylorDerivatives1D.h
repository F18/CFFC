/* TaylorDerivatives1D.h: Header file defining the specialization of 'TaylorDerivatives' template in 1D. 
   Obs: The factorial coefficients are incorporated in the values of the derivatives. */

#ifndef _TAYLORDERIVATIVES_1D_INCLUDED
#define _TAYLORDERIVATIVES_1D_INCLUDED

#ifndef _TAYLORDERIVATIVES_TEMPLATE_INCLUDED
#include "TaylorDerivatives_Template.h"
#endif // _TAYLORDERIVATIVES_TEMPLATE_INCLUDED 

/*****************************************************************************
// Define the partial specialization of the object class DerivativeObj in 1D *
*****************************************************************************/
#ifndef __Use_Iterator__
template<class T>
class DerivativeObj<OneD,T>;

/************************************************
*     Friend Functions : DerivativeObj          *
************************************************/
template<SpaceType OneD, class T >
bool operator==(const DerivativeObj<OneD,T>& left, const DerivativeObj<OneD,T>& right);

template<SpaceType OneD, class T >
bool operator!=(const DerivativeObj<OneD,T>& left, const DerivativeObj<OneD,T>& right);

template<SpaceType OneD, class T >
std::ostream& operator<< (std::ostream& os, const DerivativeObj<OneD,T>& Obj);

/*******************************************************
 * TEMPLETIZED CLASS: DerivativeObj                    *
 ******************************************************/
template< class T>
class DerivativeObj<OneD,T>{

 private:
 int Power;    		        /* the power coefficients */
 T ValueD;         		/* the Value of the Derivative, T = solution class */

 public:
 // Constructors
 DerivativeObj(void): Power(0), ValueD() { };
 DerivativeObj(const double Val): Power(0), ValueD(Val) { };
 DerivativeObj(const std::vector<unsigned> & PPP_, const T & ValueD_){
#ifndef __No_Checking__
   require(PPP_.size() == 1,"Constructor DerivativeObj failed due to difference in dimensions.\n");
#endif
   Power = PPP_[0];
   ValueD = ValueD_;
 };
 DerivativeObj(const unsigned Power_, const T ValueD_ ){
   Power = Power_;
   ValueD = ValueD_;
 }
 // Copy constructor
 DerivativeObj( const DerivativeObj & rhs): Power(rhs.Power), ValueD(rhs.ValueD){ };
 // ~DerivativeObj(){ }; use default destructor
 
 // Member functions
 void SetPowers(const unsigned p1) {Power = p1;}	/* set p1 */
 void SetPowers(const unsigned &p1, const bool) {Power = p1;} /* set p1, passing by reference */
 void SetPowers(const std::vector<unsigned> & PPP_); /* set the powers */
 void SetValue(const T ValueD_) { ValueD = ValueD_; } /* set ValueD */
 void SetValue(const T & ValueD_, const bool ) { ValueD = ValueD_;} /* set ValueD, passing by reference */
 bool IsPowerEqualTo(const unsigned p1) const ; /* if the argument provided is equal to p1 returns true */

 // Access functions
 T & D(void) { return ValueD; } 		      /* return ValueD */
 const T & D(void) const {return ValueD; }
 const double & D(const unsigned VarPosition) const { return ValueD[VarPosition]; }
 double & D(const unsigned VarPosition){ return ValueD[VarPosition]; }
 unsigned P1(void) { return Power; }	/* returns p1 */
 unsigned P1(void) const {return Power;}

 /* Overloaded operators */
 // Assignment operator
 DerivativeObj& operator=(const DerivativeObj<OneD,T>& right);
 
 // Friend functions & operators
 friend bool operator== <OneD,T> (const DerivativeObj<OneD,T>& left,
				  const DerivativeObj<OneD,T>& right);
 friend bool operator!= <OneD,T> (const DerivativeObj<OneD,T>& left,
				  const DerivativeObj<OneD,T>& right);
 friend std::ostream& operator<< <OneD,T> (std::ostream& os, const DerivativeObj<OneD,T>& Obj);
};

// CLASS DerivativeObj
// Assignment operator
template< class T>
inline DerivativeObj<OneD,T> & DerivativeObj<OneD,T>::operator=(const DerivativeObj<OneD,T>& rhs){
  if(this == &rhs) return *this;
  Power = rhs.Power;
  ValueD = rhs.ValueD;
  return *this;
}

// SetPowers()
template< class T>
inline void DerivativeObj<OneD,T>::SetPowers(const std::vector<unsigned> & PPP_){
#ifndef __No_Checking__
  require(PPP_.size() >= 0, "DerivativeObj::SetPowers(vector<unsigned>) failed. Not enough elements.\n");
#endif
  Power = PPP_[0];
}

// IsPowerEqualTo()
template< class T>
inline bool DerivativeObj<OneD,T>::IsPowerEqualTo(const unsigned p1) const {
  return (Power==p1);
}

// Friend functions
template< class T> inline
std::ostream& operator<<(std::ostream& os, const DerivativeObj<OneD,T>& Obj){
  return os << "Power=" << Obj.P1() << "\tValueD=" << Obj.D() << std::endl;
}

// operator ==
template< class T> inline
bool operator==(const DerivativeObj<OneD,T>& left,
		const DerivativeObj<OneD,T>& right){
  return (left.P1()==right.P1())&&(left.D()==right.D());
}

// operator !=
template< class T> inline
bool operator!=(const DerivativeObj<OneD,T>& left,
		const DerivativeObj<OneD,T>& right){
  return !(left == right);
}

/*******************************************************************
 * TEMPLATIZED CLASS: TaylorDerivatives 1D Container Declaration   *
 ******************************************************************/
template<class T>
class TaylorDerivativesContainer<OneD,T>;

/*******************************************************************
 * TaylorDerivativesContainer: Friend Function Declaration
 ******************************************************************/

template<SpaceType OneD, class T>
bool operator==(const TaylorDerivativesContainer<OneD,T>& left,
		const TaylorDerivativesContainer<OneD,T>& right);

template<SpaceType OneD, class T>
bool operator!= (const TaylorDerivativesContainer<OneD,T>& left,
		 const TaylorDerivativesContainer<OneD,T>& right);

template<SpaceType OneD, class T>
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<OneD,T>& Obj);

/*******************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer_1D    *
 ******************************************************/
template<class T>
class TaylorDerivativesContainer<OneD,T>{
 public:
  typedef DerivativeObj<OneD,T> Derivative;

 private:
  // Variables
  Derivative* DContainer; /* Container of derivative objects */
  int ContainerSize;
  T phi;			/* Limiter value */

  public:
  // Default Constructor
  TaylorDerivativesContainer(void): DContainer(NULL), ContainerSize(0), phi(1.0){ };
  // Overloaded constructor
  TaylorDerivativesContainer(const unsigned OrderOfReconstruction): DContainer(NULL), ContainerSize(0), phi(1.0){
    GenerateContainer(OrderOfReconstruction);
  }
  // Copy constructor
  TaylorDerivativesContainer(const TaylorDerivativesContainer<OneD,T> & rhs);
  // Destructor
  ~TaylorDerivativesContainer(void){ free_memory();}
  // Size()
  unsigned size() const { return ContainerSize;}
  unsigned FirstElem() const { return 0; }
  int LastElem() const { return ContainerSize-1;}
  int RecOrder() const { return ContainerSize-1;}
  // Allocate Memory
  void allocate(const unsigned NumberOfObjects);
  // GenerateContainer()
  void GenerateContainer(const unsigned OrderOfReconstruction); /* Allocate memory for the derivatives and 
								   generate the power combinations */
  // free_memory()
  void free_memory(); 		/* Delete the objects of the container */

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX);

  /* Field access */
  const T & Limiter(void) const {return phi;}
  T & Limiter(void) {return phi;}

  /* Overloaded Operators */
  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) {return DContainer[position];}
  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) const {return DContainer[position];}
  T & operator()(const unsigned & p1);
  T & operator()(const unsigned & p1) const;
  double & operator()(const unsigned p1, const unsigned VarPosition);
  double & operator()(const unsigned p1, const unsigned VarPosition) const;

  /* Friend functions */
  friend bool operator== <OneD,T> (const TaylorDerivativesContainer<OneD,T>& left,
			      const TaylorDerivativesContainer<OneD,T>& right);

  friend bool operator!= <OneD,T> (const TaylorDerivativesContainer<OneD,T>& left,
			      const TaylorDerivativesContainer<OneD,T>& right);

  friend std::ostream & operator<< <OneD,T> (std::ostream & os, const TaylorDerivativesContainer<OneD,T>& Obj);

    // Assignment operator;
  TaylorDerivativesContainer<OneD,T> & operator=(const TaylorDerivativesContainer<OneD,T> & rhs);

};

/*******************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer       *
 * Implementation of the Member Functions              *
 ******************************************************/

 /* Copy constructor  */
template<class T> inline
TaylorDerivativesContainer<OneD,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<OneD,T> & rhs)
{

  // allocate memory for the new container
  allocate(rhs.size());

  // copy the value from the LHS
  for (int i=0; i<=LastElem() ; ++i){
    DContainer[i] = rhs.DContainer[i];
  }
}

/* Allocate memory for the derivatives */
template<class T> inline
void TaylorDerivativesContainer<OneD,T>::allocate(const unsigned NumberOfObjects)
{
  ContainerSize = NumberOfObjects;

  if (NumberOfObjects > 0){
    /* create memory */
    DContainer = new Derivative[ContainerSize];
    
    /* initialize to ZERO */
    Derivative ZERO_Derivative(0.0);
    for(int i=0; i<ContainerSize; ++i){
      DContainer[i] = ZERO_Derivative;
    }
  } else {
    DContainer = NULL;
  }
}

/* Allocate memory for the derivatives and power coefficients */
template<class T> inline
void TaylorDerivativesContainer<OneD,T>::GenerateContainer(const unsigned OrderOfReconstruction)
{

  /* free the memory if there is memory allocated to the pointer */
  free_memory();

  /* allocate memory for the derivatives */
  allocate(OrderOfReconstruction+1); // NumberOfDerivatives = OrderOfReconstruction+1

  /* Set powers */
  for (int p1=0; p1<=OrderOfReconstruction; ++p1){
    DContainer[p1].SetPowers(p1);
  }
}

/* free_memory() */
template<class T> inline
void TaylorDerivativesContainer<OneD,T>::free_memory() {

  // Check if DContainer is empty
  if (DContainer != NULL){
    delete [] DContainer;
    DContainer = NULL;
  }
}

/* Assignment operator = */
template<class T> inline
TaylorDerivativesContainer<OneD,T> & 
TaylorDerivativesContainer<OneD,T>::operator=(const TaylorDerivativesContainer<OneD,T> & rhs){

  // !!! If the LHS container already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & rhs) return *this;

  // Check if there is memory allocated
  free_memory();

  // allocate memory
  allocate(rhs.size());

  // copy the value from the LHS
  for (int i=0; i<=LastElem() ; ++i){
    DContainer[i] = rhs.DContainer[i];
  }

  return *this;
}

/* Access data */
// operator(unsigned) -> returns the derivative for which p1_ = p1
template<class T> inline
T & TaylorDerivativesContainer<OneD,T>::operator()(const unsigned & p1){

  return DContainer[p1].D();
}

// operator(unsigned) -> returns the derivative for which p1_ = p1
template<class T> inline 
T & TaylorDerivativesContainer<OneD,T>::operator()(const unsigned & p1) const{

  return DContainer[p1].D();
}

template<class T> inline
double & TaylorDerivativesContainer<OneD,T>::operator()(const unsigned p1, const unsigned VarPosition){
  return DContainer[p1].D(VarPosition);
}

template<class T> inline
double & TaylorDerivativesContainer<OneD,T>::operator()(const unsigned p1, const unsigned VarPosition) const{
  return DContainer[p1].D(VarPosition);
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance DeltaX
template<class T> inline
T TaylorDerivativesContainer<OneD,T>::ComputeSolutionFor(const double DeltaX){

  // initialize the solution state
  T Solution(0.0);
  double DeltaXtoPower(1.0);

  // compute the solution using the derivatives
  for (int p1=1; p1<=LastElem(); ++p1){
    /* compute DeltaX to power "p1" */
    DeltaXtoPower *= DeltaX;
    /* add the contribution of the p1 derivative */
    Solution += DeltaXtoPower * DContainer[p1].D();
  }

  return DContainer[0].D() + phi*Solution;
}

// Friend functions
template<class T> inline
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<OneD,T>& Obj){

  for(int i=0; i<=Obj.LastElem(); ++i)
    os << Obj(i,true,true,true);
  return os << std::endl;
}

// operator ==
template<class T> inline
bool operator==(const TaylorDerivativesContainer<OneD,T>& left,
		const TaylorDerivativesContainer<OneD,T>& right){

  bool answer = true;
  if (left.size() != right.size() )
    return !answer;
  for (int i=0; i<=right.LastElem(); ++i){
    answer = answer && ( left(i,true,true,true) == right(i,true,true,true) );
  }
  return answer;
}

// operator !=
template<class T>
inline bool operator!=(const TaylorDerivativesContainer<OneD,T>& left,
		       const TaylorDerivativesContainer<OneD,T>& right){

  return !(left == right);
}

#endif // __Use_Iterator__

// Specialization for "class T == double"
template<> inline
const double & DerivativeObj<OneD,double>::D(const unsigned ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<OneD,double>::D(const unsigned ) {
  return ValueD;
}

#endif
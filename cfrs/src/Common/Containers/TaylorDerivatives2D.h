/* TaylorDerivatives1D.h: Header file defining the specialization of 'TaylorDerivatives' template in 1D. 
   Obs: The factorial coefficients are incorporated in the values of the derivatives. */

#ifndef _TAYLORDERIVATIVES_2D_INCLUDED
#define _TAYLORDERIVATIVES_2D_INCLUDED

#ifndef _TAYLORDERIVATIVES_TEMPLATE_INCLUDED
#include "TaylorDerivatives_Template.h"
#endif // _TAYLORDERIVATIVES_TEMPLATE_INCLUDED 

/*****************************************************************************
// Define the partial specialization of the object class DerivativeObj in 2D *
*****************************************************************************/
#ifndef __Use_Iterator__
template<class T>
class DerivativeObj<TwoD,T>;

/************************************************
*     Friend Functions : DerivativeObj          *
************************************************/
template<SpaceType TwoD, class T >
bool operator==(const DerivativeObj<TwoD,T>& left, const DerivativeObj<TwoD,T>& right);

template<SpaceType TwoD, class T >
bool operator!=(const DerivativeObj<TwoD,T>& left, const DerivativeObj<TwoD,T>& right);

template<SpaceType TwoD, class T >
std::ostream& operator<< (std::ostream& os, const DerivativeObj<TwoD,T>& Obj);

template<SpaceType TwoD, class T >
std::istream& operator>> (std::istream& os, DerivativeObj<TwoD,T>& Obj);

/*******************************************************
 * CLASS Template: DerivativeObj                    *
 ******************************************************/
template< class T>
class DerivativeObj<TwoD,T>{

 private:
 unsigned Power1;  		        /* the power coefficients */
 unsigned Power2;
 T ValueD;         		/* the Value of the Derivative, T = solution class */

 public:
 // Constructors
 DerivativeObj(void): Power1(0), Power2(0), ValueD() { };
 DerivativeObj(const double Val): Power1(0), Power2(0), ValueD(Val) { };
 DerivativeObj(const std::vector<int> & PPP_, const T & ValueD_);
 DerivativeObj(const unsigned p1_, const unsigned p2_, const T ValueD_): Power1(p1_),Power2(p2_),ValueD(ValueD_) { };
 // Copy constructor
 DerivativeObj( const DerivativeObj & rhs): Power1(rhs.Power1), Power2(rhs.Power2), ValueD(rhs.ValueD){ };
 ~DerivativeObj(){ };
 
 // Member functions
 void SetPowers(const unsigned p1, const unsigned p2) {Power1 = p1; Power2 = p2; }	/* set p1 & p2 */
 /* set p1 & p2 , passing by reference */
 void SetPowers(const unsigned &p1, const unsigned &p2, const bool) {Power1 = p1; Power2 = p2;}
 void SetPowers(const std::vector<unsigned> & PPP_); /* set the powers */
 void SetValue(const T ValueD_) { ValueD = ValueD_; } /* set ValueD */
 void SetValue(const T & ValueD_, const bool ) { ValueD = ValueD_;} /* set ValueD, passing by reference */
 bool IsPowerEqualTo(const unsigned p1, const unsigned p2) const ; /* if the argument provided is equal to p1 & p2 returns true */

 void Read_Derivative(istream &In_File);

 // Access functions
 T & D(void) { return ValueD; } 		      /* return ValueD */
 const T & D(void) const {return ValueD; }
 const double & D(const unsigned VarPosition) const { return ValueD[VarPosition]; }
 double & D(const unsigned VarPosition){ return ValueD[VarPosition]; }
 unsigned P1(void) { return Power1; }	/* returns p1 */
 unsigned P1(void) const {return Power1;}
 unsigned P2(void) { return Power2; }	/* returns p2 */
 unsigned P2(void) const {return Power2;}

 /* Overloaded operators */
 // Assignment operator
 DerivativeObj& operator=(const DerivativeObj<TwoD,T>& right);
 
 // Friend functions & operators
 friend bool operator== <TwoD,T> (const DerivativeObj<TwoD,T>& left,
				  const DerivativeObj<TwoD,T>& right);
 friend bool operator!= <TwoD,T> (const DerivativeObj<TwoD,T>& left,
				  const DerivativeObj<TwoD,T>& right);
 friend std::ostream& operator<< <TwoD,T> (std::ostream& os, const DerivativeObj<TwoD,T>& Obj);

 friend std::istream& operator>> <TwoD,T> (std::istream& os, DerivativeObj<TwoD,T>& Obj);
};

// CLASS DerivativeObj

// Constructor
template <class T> inline
DerivativeObj<TwoD,T>::DerivativeObj(const std::vector<int> & PPP_, const T & ValueD_){
#ifndef __No_Checking__
   require(PPP_.size() == 2,"Constructor DerivativeObj failed due to difference in dimensions.\n");
#endif
   Power1 = PPP_[0];
   Power2 = PPP_[1];
   ValueD = ValueD_;
}

// Assignment operator
template< class T>
inline DerivativeObj<TwoD,T> & DerivativeObj<TwoD,T>::operator=(const DerivativeObj<TwoD,T>& rhs){
  if(this == &rhs) return *this;
  Power1 = rhs.Power1;
  Power2 = rhs.Power2;
  ValueD = rhs.ValueD;
  return *this;
}

// SetPowers()
template< class T>
inline void DerivativeObj<TwoD,T>::SetPowers(const std::vector<unsigned> & PPP_){
#ifndef __No_Checking__
  require(PPP_.size() >= 1, "DerivativeObj::SetPowers(vector<unsigned>) failed. Not enough elements.\n");
#endif
  Power1 = PPP_[0];
  Power2 = PPP_[1];
}

// IsPowerEqualTo()
template< class T>
inline bool DerivativeObj<TwoD,T>::IsPowerEqualTo(const unsigned p1, const unsigned p2) const {
  return ( Power1==p1 && Power2==p2 );
}

// Friend functions
template< class T> inline
std::ostream& operator<<(std::ostream& os, const DerivativeObj<TwoD,T>& Obj){
  return os << "powers= " << Obj.P1() << ", " << Obj.P2() << ", ValueD=" << Obj.D() << std::endl;
}

template< class T> inline
void DerivativeObj<TwoD,T>::Read_Derivative(istream &In_File){

  In_File.setf(ios::skipws);
  In_File >> Power1 >> Power2 >> ValueD; 
  In_File.unsetf(ios::skipws);
}

template< class T>
std::istream& operator>>(std::istream& os, DerivativeObj<TwoD,T>& Obj){

  Obj.Read_Derivative(os);

  return os;
}

// operator ==
template< class T> inline
bool operator==(const DerivativeObj<TwoD,T>& left,
		const DerivativeObj<TwoD,T>& right){
  return (left.P1()==right.P1())&&(left.P2()==right.P2())&&(left.D()==right.D());
}

// operator !=
template< class T> inline
bool operator!=(const DerivativeObj<TwoD,T>& left,
		const DerivativeObj<TwoD,T>& right){
  return !(left == right);
}

/*******************************************************************
 * CLASS Templete: TaylorDerivatives 2D Container Declaration   *
 ******************************************************************/
template<class T>
class TaylorDerivativesContainer<TwoD,T>;

/*******************************************************************
 * TaylorDerivativesContainer: Friend Function Declaration
 ******************************************************************/

template<SpaceType TwoD, class T>
bool operator==(const TaylorDerivativesContainer<TwoD,T>& left,
		const TaylorDerivativesContainer<TwoD,T>& right);

template<SpaceType TwoD, class T>
bool operator!= (const TaylorDerivativesContainer<TwoD,T>& left,
		 const TaylorDerivativesContainer<TwoD,T>& right);

template<SpaceType TwoD, class T>
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj);

template<SpaceType TwoD, class T>
std::istream & operator>> (std::istream & os, TaylorDerivativesContainer<TwoD,T>& Obj);

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer_2D    *
 ******************************************************/
template<class T>
class TaylorDerivativesContainer<TwoD,T>{
 public:
  typedef DerivativeObj<TwoD,T> Derivative;

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
  TaylorDerivativesContainer(const TaylorDerivativesContainer<TwoD,T> & rhs);
  // Destructor
  ~TaylorDerivativesContainer(void){ free_memory();}
  // Assignment operator;
  TaylorDerivativesContainer<TwoD,T> & operator=(const TaylorDerivativesContainer<TwoD,T> & rhs);

  // Memory Management
  void allocate(const unsigned NumberOfObjects);
  // free_memory()
  void free_memory(); 		/* Delete the objects of the container */
  // GenerateContainer()
  void GenerateContainer(const unsigned OrderOfReconstruction); /* Allocate memory for the derivatives and 
								   generate the power combinations */

  // Characteristic parameters
  unsigned size() const { return ContainerSize;}
  unsigned FirstElem() const { return 0; }
  int LastElem() const { return ContainerSize-1;}
  unsigned RecOrder() const { return (unsigned)(sqrt(1.0+8.0*ContainerSize)-3.0)/2.0;}
  // IndexOrder(unsigned,unsigned) -> determines the position of the pointer to the element having the powers (p1,p2)
  unsigned IndexOrder(const unsigned p1, const unsigned p2);

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX, const double DeltaY);

  /* Field access */
  const T & Limiter(void) const {return phi;}
  T & Limiter(void) {return phi;}

  /* Overloaded Operators */
  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) {return DContainer[position];}
  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) const {return DContainer[position];}
  T & operator()(const unsigned p1, const unsigned p2);
  T & operator()(const unsigned p1, const unsigned p2) const;
  double & operator()(const unsigned p1, const unsigned p2, const unsigned VarPosition);
  double & operator()(const unsigned p1, const unsigned p2, const unsigned VarPosition) const;

  /* Friend functions */
  friend bool operator== <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
				   const TaylorDerivativesContainer<TwoD,T>& right);

  friend bool operator!= <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
				   const TaylorDerivativesContainer<TwoD,T>& right);

  friend std::ostream & operator<< <TwoD,T> (std::ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj);

  friend std::istream & operator>> <TwoD,T> (std::istream & os, TaylorDerivativesContainer<TwoD,T>& Obj);

};

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer          *
 * Implementation of the Member Functions              *
 ******************************************************/

 /* Copy constructor  */
template<class T> inline
TaylorDerivativesContainer<TwoD,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<TwoD,T> & rhs)
{

  // allocate memory for the new container
  allocate(rhs.size());

  // copy the value from the LHS
  for (int i=0; i<=LastElem(); ++i){
    DContainer[i] = rhs.DContainer[i];
  }
}

/* Allocate memory for the derivatives */
template<class T> inline
void TaylorDerivativesContainer<TwoD,T>::allocate(const unsigned NumberOfObjects)
{

  ContainerSize = NumberOfObjects;

  if (ContainerSize > 0){
    /* create memory */    
    DContainer = new Derivative[ContainerSize];
    
    /* initialize state to ZERO */
    Derivative ZERO_Derivative(0.0);
    for(int i=0; i<=LastElem(); ++i){
      DContainer[i] = ZERO_Derivative;
    }
  } else {
    DContainer = NULL;
  }
}

/* Allocate memory for the derivatives and power coefficients */
template<class T> inline
void TaylorDerivativesContainer<TwoD,T>::GenerateContainer(const unsigned OrderOfReconstruction)
{

  int NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)/2;

  /* free the memory if there is memory allocated to the pointer */
  free_memory();

  /* allocate memory for the derivatives */
  allocate(NumberOfDerivatives); // 

  int p1,p2,Position;

  /* Set powers */
  for (p1=0,Position=0; p1<=OrderOfReconstruction; ++p1)
    for (p2=0; p2<=OrderOfReconstruction-p1; ++p2, ++Position){
      DContainer[Position].SetPowers(p1,p2);
    }
}

/* free_memory() */
template<class T> inline
void TaylorDerivativesContainer<TwoD,T>::free_memory() {

  // Check if DContainer is empty
  if (DContainer != NULL){
    delete [] DContainer;
    DContainer = NULL;
  }
}

/* Assignment operator = */
template<class T> inline
TaylorDerivativesContainer<TwoD,T> & 
TaylorDerivativesContainer<TwoD,T>::operator=(const TaylorDerivativesContainer<TwoD,T> & rhs){

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
/*IndexOrder(unsigned,unsigned) -> determines the position of the element having the power combination (p1,p2) */
template<class T> inline
  unsigned TaylorDerivativesContainer<TwoD,T>::IndexOrder(const unsigned p1, const unsigned p2)
{

  /* Obs. DContainer[LastElem()].P1() is equal to the order of the reconstruction */
  unsigned sum = 0;

  for (unsigned n=2; n<=p1; ++n){
    sum += (n - 1);
  }

  return (DContainer[LastElem()].P1() + 1)*p1  - sum + p2;
}


// operator(unsigned) -> returns the derivative for which Power1 = p1 & Power2 = p2
template<class T> inline
T & TaylorDerivativesContainer<TwoD,T>::operator()(const unsigned p1, const unsigned p2){

  return DContainer[IndexOrder(p1,p2)].D();
}

// operator(unsigned) -> returns the derivative for which p1_ = p1
template<class T> inline 
T & TaylorDerivativesContainer<TwoD,T>::operator()(const unsigned p1, const unsigned p2) const{

  return DContainer[IndexOrder(p1,p2)].D();
}

template<class T> inline
double & TaylorDerivativesContainer<TwoD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned VarPosition){
  return DContainer[IndexOrder(p1,p2)].D(VarPosition);
}

template<class T> inline
double & TaylorDerivativesContainer<TwoD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned VarPosition) const{
  return DContainer[IndexOrder(p1,p2)].D(VarPosition);
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY)
template<class T> inline
T TaylorDerivativesContainer<TwoD,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY){

  // initialize the solution state
  T Solution(0.0);
  T OneT(1.0);

  int OrderOfReconstruction = DContainer[LastElem()].P1();
  int p1,p2,Position;

  double DeltaXtoPower(1.0), DeltaYtoPower;

  for (p1=0,Position=0; p1<=OrderOfReconstruction; ++p1){
    /* Reinitialize DeltaYtoPower */
    DeltaYtoPower = 1.0;
    for (p2=0; p2<=OrderOfReconstruction-p1; ++p2, ++Position){
      /* Update solution */
      Solution += DeltaXtoPower*DeltaYtoPower*DContainer[Position].D();

      /* Update DeltaYtoPower */
      DeltaYtoPower *= DeltaY;
    }

    /* Update DeltaXtoPower */
    DeltaXtoPower *= DeltaX;
  }

  return phi*Solution + (OneT-phi)*DContainer[0].D();
}

// Friend functions
template<class T> inline
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj){

  for(int i=0; i<=Obj.LastElem(); ++i)
    os << Obj(i,true,true,true);
  return os << std::endl;
}

template<class T> inline
std::istream & operator>> (std::istream & os, TaylorDerivativesContainer<TwoD,T>& Obj){

  for(int i=0; i<=Obj.LastElem(); ++i){
    os >> Obj(i,true,true,true);
  }
  return os;
}

// operator ==
template<class T> inline
bool operator==(const TaylorDerivativesContainer<TwoD,T>& left,
		const TaylorDerivativesContainer<TwoD,T>& right){

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
inline bool operator!=(const TaylorDerivativesContainer<TwoD,T>& left,
		       const TaylorDerivativesContainer<TwoD,T>& right){

  return !(left == right);
}

#endif // __Use_Iterator__

// Specialization for "class T == double"
template<> inline
const double & DerivativeObj<TwoD,double>::D(const unsigned ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<TwoD,double>::D(const unsigned ) {
  return ValueD;
}

#endif

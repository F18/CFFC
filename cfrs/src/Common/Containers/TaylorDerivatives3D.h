/* TaylorDerivatives3D.h: Header file defining the specialization of 'TaylorDerivatives' template in 3D. 
   Obs: The factorial coefficients are incorporated in the values of the derivatives. */

#ifndef _TAYLORDERIVATIVES_3D_INCLUDED
#define _TAYLORDERIVATIVES_3D_INCLUDED

#ifndef _TAYLORDERIVATIVES_TEMPLATE_INCLUDED
#include "TaylorDerivatives_Template.h"
#endif // _TAYLORDERIVATIVES_TEMPLATE_INCLUDED 

/*****************************************************************************
// Define the partial specialization of the object class DerivativeObj in 3D *
*****************************************************************************/
#ifndef __Use_Iterator__
template<class T>
class DerivativeObj<ThreeD,T>;

/************************************************
*     Friend Functions : DerivativeObj          *
************************************************/
template<SpaceType ThreeD, class T >
bool operator==(const DerivativeObj<ThreeD,T>& left, const DerivativeObj<ThreeD,T>& right);

template<SpaceType ThreeD, class T >
bool operator!=(const DerivativeObj<ThreeD,T>& left, const DerivativeObj<ThreeD,T>& right);

template<SpaceType ThreeD, class T >
std::ostream& operator<< (std::ostream& os, const DerivativeObj<ThreeD,T>& Obj);

template<SpaceType ThreeD, class T >
std::istream& operator>> (std::istream& os, DerivativeObj<ThreeD,T>& Obj);

/*******************************************************
 * CLASS Template: DerivativeObj                    *
 ******************************************************/
template< class T>
class DerivativeObj<ThreeD,T>{
  
private:
  unsigned Power1;  		        /* the power coefficients */
  unsigned Power2;
  unsigned Power3;
  T ValueD;         		/* the Value of the Derivative, T = solution class */
  
public:
  // Constructors
  DerivativeObj(void): Power1(0), Power2(0), Power3(0), ValueD() { };
  DerivativeObj(const double Val): Power1(0), Power2(0), Power3(0), ValueD(Val) { };
  DerivativeObj(const std::vector<int> & PPP_, const T & ValueD_);
  DerivativeObj(const unsigned p1_, const unsigned p2_, const unsigned p3_, const T ValueD_): Power1(p1_),Power2(p2_),Power3(p3_),ValueD(ValueD_) { };
  // Copy constructor
  DerivativeObj( const DerivativeObj & rhs): Power1(rhs.Power1), Power2(rhs.Power2), Power3(rhs,Power3),ValueD(rhs.ValueD){ };
  ~DerivativeObj(){ };
 
  // Member functions
  void SetPowers(const unsigned p1, const unsigned p2, const unsigned p3) {Power1 = p1; Power2 = p2; Power3 = p3; }	/* set p1, p2, & p3 */
  /* set p1, p2, & p3 , passing by reference */
  void SetPowers(const unsigned &p1, const unsigned &p2, const unsigned &p3, const bool) {Power1 = p1; Power2 = p2; Power3=p3;}
  void SetPowers(const std::vector<unsigned> & PPP_); /* set the powers */
  void SetValue(const T ValueD_) { ValueD = ValueD_; } /* set ValueD */
  void SetValue(const T & ValueD_, const bool ) { ValueD = ValueD_;} /* set ValueD, passing by reference */
  bool IsPowerEqualTo(const unsigned p1, const unsigned p2, const unsigned p3) const ; /* if the argument provided is equal to p1, p2, & p3 returns true */

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
  unsigned P3(void) { return Power3; }	/* returns p3 */
  unsigned P3(void) const {return Power3;}

  /* Overloaded operators */
  // Assignment operator
  DerivativeObj& operator=(const DerivativeObj<ThreeD,T>& right);
 
  // Friend functions & operators
  friend bool operator== <ThreeD,T> (const DerivativeObj<ThreeD,T>& left,
				     const DerivativeObj<ThreeD,T>& right);
  friend bool operator!= <ThreeD,T> (const DerivativeObj<ThreeD,T>& left,
				     const DerivativeObj<ThreeD,T>& right);
  friend std::ostream& operator<< <ThreeD,T> (std::ostream& os, const DerivativeObj<ThreeD,T>& Obj);
  
  friend std::istream& operator>> <ThreeD,T> (std::istream& os, DerivativeObj<ThreeD,T>& Obj);
};

// CLASS DerivativeObj

// Constructor
template <class T> inline
DerivativeObj<ThreeD,T>::DerivativeObj(const std::vector<int> & PPP_, const T & ValueD_){
#ifndef __No_Checking__
   require(PPP_.size() == 3,"Constructor DerivativeObj failed due to difference in dimensions.\n");
#endif
   Power1 = PPP_[0];
   Power2 = PPP_[1];
   Power3 = PPP_[2];
   ValueD = ValueD_;
}

// Assignment operator
template< class T>
inline DerivativeObj<ThreeD,T> & DerivativeObj<ThreeD,T>::operator=(const DerivativeObj<ThreeD,T>& rhs){
  if(this == &rhs) return *this;
  Power1 = rhs.Power1;
  Power2 = rhs.Power2;
  Power3 = rhs.Power3;
  ValueD = rhs.ValueD;
  return *this;
}

// SetPowers()
template< class T>
inline void DerivativeObj<ThreeD,T>::SetPowers(const std::vector<unsigned> & PPP_){
#ifndef __No_Checking__
  require(PPP_.size() >= 3, "DerivativeObj::SetPowers(vector<unsigned>) failed. Not enough elements.\n");
#endif
  Power1 = PPP_[0];
  Power2 = PPP_[1];
  Power3 = PPP_[2];
}

// IsPowerEqualTo()
template< class T>
inline bool DerivativeObj<ThreeD,T>::IsPowerEqualTo(const unsigned p1, const unsigned p2, const unsigned p3) const {
  return ( Power1==p1 && Power2==p2 && Power3==p3 );
}

// Friend functions
template< class T> inline
std::ostream& operator<<(std::ostream& os, const DerivativeObj<ThreeD,T>& Obj){
  return os << "powers= " << Obj.P1() << ", " << Obj.P2() << ", " << Obj.P3() << ", ValueD= " << Obj.D() << std::endl;
}

template< class T> inline
void DerivativeObj<ThreeD,T>::Read_Derivative(istream &In_File){

  In_File.setf(ios::skipws);
  In_File >> Power1 >> Power2 >> Power3 >> ValueD; 
  In_File.unsetf(ios::skipws);
}

template< class T>
std::istream& operator>>(std::istream& os, DerivativeObj<ThreeD,T>& Obj){

  Obj.Read_Derivative(os);

  return os;
}

// operator ==
template< class T> inline
bool operator==(const DerivativeObj<ThreeD,T>& left,
		const DerivativeObj<ThreeD,T>& right){
  return ( (left.P1()==right.P1()) && (left.P2()==right.P2()) && (left.P3()==right.P3()) && (left.D()==right.D()) );
}

// operator !=
template< class T> inline
bool operator!=(const DerivativeObj<ThreeD,T>& left,
		const DerivativeObj<ThreeD,T>& right){
  return !(left == right);
}

/*******************************************************************
 * CLASS Template: TaylorDerivatives 3D Container Declaration   *
 ******************************************************************/
template<class T>
class TaylorDerivativesContainer<ThreeD,T>;

/*******************************************************************
 * TaylorDerivativesContainer: Friend Function Declaration
 ******************************************************************/

template<SpaceType ThreeD, class T>
bool operator==(const TaylorDerivativesContainer<ThreeD,T>& left,
		const TaylorDerivativesContainer<ThreeD,T>& right);

template<SpaceType ThreeD, class T>
bool operator!= (const TaylorDerivativesContainer<ThreeD,T>& left,
		 const TaylorDerivativesContainer<ThreeD,T>& right);

template<SpaceType ThreeD, class T>
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj);

template<SpaceType ThreeD, class T>
std::istream & operator>> (std::istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj);

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer_3D    *
 ******************************************************/
template<class T>
class TaylorDerivativesContainer<ThreeD,T>{
 public:
  typedef DerivativeObj<ThreeD,T> Derivative;

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
  TaylorDerivativesContainer(const TaylorDerivativesContainer<ThreeD,T> & rhs);
  // Destructor
  ~TaylorDerivativesContainer(void){ free_memory();}
  // Assignment operator;
  TaylorDerivativesContainer<ThreeD,T> & operator=(const TaylorDerivativesContainer<ThreeD,T> & rhs);

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
  // IndexOrder(unsigned,unsigned,unsigned) -> determines the position of the pointer to the element having the powers (p1,p2,p3)
  unsigned IndexOrder(const unsigned p1, const unsigned p2, const unsigned p3);

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ);

  /* Field access */
  const T & Limiter(void) const {return phi;}
  T & Limiter(void) {return phi;}

  /* Overloaded Operators */
  //  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) {return DContainer[position];}
  //  Derivative & operator()(const unsigned & position, const bool, const bool, const bool) const {return DContainer[position];}
  Derivative & operator()(const unsigned & position) {return DContainer[position];}
  Derivative & operator()(const unsigned & position) const {return DContainer[position];}
  T & operator()(const unsigned p1, const unsigned p2, const unsigned p3);
  T & operator()(const unsigned p1, const unsigned p2, const unsigned p3) const;
  double  & operator()(const unsigned p1, const unsigned p2, const unsigned p3, const unsigned VarPosition);
  double  & operator()(const unsigned p1, const unsigned p2, const unsigned p3, const unsigned VarPosition) const;

  /* Friend functions */
  friend bool operator== <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
				   const TaylorDerivativesContainer<ThreeD,T>& right);

  friend bool operator!= <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
				   const TaylorDerivativesContainer<ThreeD,T>& right);

  friend std::ostream & operator<< <ThreeD,T> (std::ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj);

  friend std::istream & operator>> <ThreeD,T> (std::istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj);

};

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer          *
 * Implementation of the Member Functions              *
 ******************************************************/

 /* Copy constructor  */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<ThreeD,T> & rhs)
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
void TaylorDerivativesContainer<ThreeD,T>::allocate(const unsigned NumberOfObjects)
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
void TaylorDerivativesContainer<ThreeD,T>::GenerateContainer(const unsigned OrderOfReconstruction)
{

  // int NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)/2
  int NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)*(OrderOfReconstruction+3)/6;

  /* free the memory if there is memory allocated to the pointer */
  free_memory();

  /* allocate memory for the derivatives */
  allocate(NumberOfDerivatives); // 

  int p1,p2,p3,Position;

  /* Set powers */
  for (p1=0,Position=0; p1<=OrderOfReconstruction; ++p1)
    for (p2=0; p2<=OrderOfReconstruction-p1; ++p2)
      for (p3=0; p3<=OrderOfReconstruction-p1-p2; ++p3, ++Position){
	DContainer[Position].SetPowers(p1,p2,p3);
    }
}

/* free_memory() */
template<class T> inline
void TaylorDerivativesContainer<ThreeD,T>::free_memory() {

  // Check if DContainer is empty
  if (DContainer != NULL){
    delete [] DContainer;
    DContainer = NULL;
  }
}

/* Assignment operator = */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T> & 
TaylorDerivativesContainer<ThreeD,T>::operator=(const TaylorDerivativesContainer<ThreeD,T> & rhs){

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
/*IndexOrder(unsigned,unsigned,unsigned) -> determines the position of the element having the power combination (p1,p2,p3) */
template<class T> inline
unsigned TaylorDerivativesContainer<ThreeD,T>::IndexOrder(const unsigned p1, const unsigned p2, const unsigned p3)
{
  /* Observe that DContainer[LastElem()].P1() is equal to the order of the reconstruction */
  int OrderOfReconstruction = DContainer[LastElem()].P1();
  int sum_series[20];
  unsigned sum_p1 = 0;
  unsigned sum_p2 = 0;
  unsigned n;

  // sum_series is a vector containing the series 1,3,6,10,15,21,28,...
  sum_series[0] = 1;
  for (n=1; n<20; ++n){
    sum_series[n] = sum_series[n-1] + (n+1);
  }

  sum_p1 = 0;
  for (n=1; n<=p1; ++n){
    sum_p1 += sum_series[OrderOfReconstruction-n+1];
  }

  sum_p2 = 0;
  for ( n=2; n<=p2; ++n){
    sum_p2 += (n-1);
  }
  sum_p2 += p1*p2;

  return sum_p1 + (OrderOfReconstruction+1)*p2 - sum_p2 + p3;
}


// operator(unsigned) -> returns the derivative for which Power1 = p1 & Power2 = p2 & Power3 = p3
template<class T> inline
T & TaylorDerivativesContainer<ThreeD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned p3){

  return DContainer[IndexOrder(p1,p2,p3)].D();
}

// operator(unsigned) -> returns the derivative for which p1_ = p1, p2_ = p2, p3_ = p3
template<class T> inline 
T & TaylorDerivativesContainer<ThreeD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned p3) const{

  return DContainer[IndexOrder(p1,p2,p3)].D();
}

template<class T> inline
double & TaylorDerivativesContainer<ThreeD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned p3, const unsigned VarPosition){
  return DContainer[IndexOrder(p1,p2,p3)].D(VarPosition);
}

template<class T> inline
double & TaylorDerivativesContainer<ThreeD,T>::operator()(const unsigned p1, const unsigned p2, const unsigned p3, const unsigned VarPosition) const{
  return DContainer[IndexOrder(p1,p2,p3)].D(VarPosition);
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<class T> inline
T TaylorDerivativesContainer<ThreeD,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ){

  // initialize the solution state
  T Solution(0.0);
  T OneT(1.0);

  int OrderOfReconstruction = DContainer[LastElem()].P1();
  int p1(0),p2(0),p3(0),Position(0);

  double DeltaXtoPower(1.0), DeltaYtoPower, DeltaZtoPower;

 for (p1=0,Position=0; p1<=OrderOfReconstruction; ++p1){
    /* Reinitialize DeltaYtoPower */
    DeltaYtoPower = 1.0;

    for (p2=0; p2<=OrderOfReconstruction-p1; ++p2){
      /* Reinitialize DeltaZtoPower */
      DeltaZtoPower = 1.0;
      
      for (p3=0; p3<=OrderOfReconstruction-p1-p2; ++p3, ++Position){

	/* Update solution */
	 Solution += DeltaXtoPower*DeltaYtoPower*DeltaZtoPower*DContainer[Position].D();
	/* Update DeltaZtoPower */
	DeltaZtoPower *= DeltaZ;
      }

      /* Update DeltaYtoPower */
      DeltaYtoPower *= DeltaY;
    }

    /* Update DeltaXtoPower */
    DeltaXtoPower *= DeltaX;
  }
 //std::cout << "ComputeSolutionFor = " << phi*Solution + (OneT-phi)*DContainer[0].D() << endl;

 return phi*Solution + (OneT-phi)*DContainer[0].D();
}

// Friend functions
template<class T> inline
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj){

  for(int i=0; i<=Obj.LastElem(); ++i)
    os << Obj(i,true,true,true);
  return os << std::endl;
}

template<class T> inline
std::istream & operator>> (std::istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj){

  for(int i=0; i<=Obj.LastElem(); ++i){
    os >> Obj(i,true,true,true);
  }
  return os;
}

// operator ==
template<class T> inline
bool operator==(const TaylorDerivativesContainer<ThreeD,T>& left,
		const TaylorDerivativesContainer<ThreeD,T>& right){

  bool answer = true;
  if (left.size() != right.size() )
    return !answer;
  for (int i=0; i<=right.LastElem(); ++i){
    answer = answer && ( left(i) == right(i) );
  }
  return answer;
}

// operator !=
template<class T>
inline bool operator!=(const TaylorDerivativesContainer<ThreeD,T>& left,
		       const TaylorDerivativesContainer<ThreeD,T>& right){

  return !(left == right);
}

#endif // __Use_Iterator__

// Specialization for "class T == double"
template<> inline
const double & DerivativeObj<ThreeD,double>::D(const unsigned ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<ThreeD,double>::D(const unsigned ) {
  return ValueD;
}

#endif

/*!\file TaylorDerivatives3D.h
  \brief Header file defining the specialization of 'TaylorDerivatives' template in 3D. 
  Obs: The factorial coefficients are incorporated in the values of the derivatives. */

#ifndef _TAYLORDERIVATIVES_3D_INCLUDED
#define _TAYLORDERIVATIVES_3D_INCLUDED

/* Include required C++ libraries. */
#include <cassert>

/* Using std namespace functions */
// None

/*****************************************************************************
// Define the partial specialization of the object class DerivativeObj in 3D *
*****************************************************************************/
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
  int Power1;  		        /* the power coefficients */
  int Power2;
  int Power3;
  T ValueD;         		/* the Value of the Derivative, T = solution class */
  
public:
  // Constructors
  DerivativeObj(void): Power1(0), Power2(0), Power3(0), ValueD() { };
  DerivativeObj(const double Val): Power1(0), Power2(0), Power3(0), ValueD(Val) { };
  DerivativeObj(const vector<int> & PPP_, const T & ValueD_);
  DerivativeObj(const int p1_, const int p2_, const int p3_, const T ValueD_): Power1(p1_),Power2(p2_),Power3(p3_),ValueD(ValueD_) { };
  // Copy constructor
  DerivativeObj( const DerivativeObj & rhs): Power1(rhs.Power1), Power2(rhs.Power2), Power3(rhs,Power3),ValueD(rhs.ValueD){ };
  ~DerivativeObj(){ };
 
  // Member functions
  void SetPowers(const int p1, const int p2, const int p3) {Power1 = p1; Power2 = p2; Power3 = p3; }	/* set p1, p2, & p3 */
  /* set p1, p2, & p3 , passing by reference */
  void SetPowers(const int &p1, const int &p2, const int &p3, const bool) {Power1 = p1; Power2 = p2; Power3=p3;}
  void SetPowers(const vector<int> & PPP_); /* set the powers */
  void SetValue(const T ValueD_) { ValueD = ValueD_; } /* set ValueD */
  void SetValue(const T & ValueD_, const bool ) { ValueD = ValueD_;} /* set ValueD, passing by reference */
  bool IsPowerEqualTo(const int p1, const int p2, const int p3) const ; /* if the argument provided is equal to p1, p2, & p3 returns true */

  void Read_Derivative(istream &In_File);

  // Access functions
  T & D(void) { return ValueD; } 		      /* return ValueD */
  const T & D(void) const {return ValueD; }
  const double & D(const int VarPosition) const { return ValueD[VarPosition]; }
  double & D(const int VarPosition){ return ValueD[VarPosition]; }
  short int P1(void) { return Power1; }	/* returns p1 */
  const short int P1(void) const {return Power1;}
  short int P2(void) { return Power2; }	/* returns p2 */
  const short int P2(void) const {return Power2;}
  short int P3(void) { return Power3; }	/* returns p3 */
  const short int P3(void) const {return Power3;}

  /* Overloaded operators */
  // Assignment operator
  DerivativeObj& operator=(const DerivativeObj<ThreeD,T>& rhs);
 
  // Friend functions & operators
  friend bool operator== <ThreeD,T> (const DerivativeObj<ThreeD,T>& left,
				     const DerivativeObj<ThreeD,T>& right);
  friend bool operator!= <ThreeD,T> (const DerivativeObj<ThreeD,T>& left,
				     const DerivativeObj<ThreeD,T>& right);
  friend ostream& operator<< <ThreeD,T> (ostream& os, const DerivativeObj<ThreeD,T>& Obj);
  
  friend istream& operator>> <ThreeD,T> (istream& os, DerivativeObj<ThreeD,T>& Obj);
};

// CLASS DerivativeObj

// Constructor
template <class T> inline
DerivativeObj<ThreeD,T>::DerivativeObj(const vector<int> & PPP_, const T & ValueD_){
   require(PPP_.size() == 3,"Constructor DerivativeObj failed due to difference in dimensions.\n");
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
inline void DerivativeObj<ThreeD,T>::SetPowers(const vector<int> & PPP_){
  require(PPP_.size() >= 3, "DerivativeObj::SetPowers(vector<int>) failed. Not enough elements.\n");
  Power1 = PPP_[0];
  Power2 = PPP_[1];
  Power3 = PPP_[2];
}

// IsPowerEqualTo()
template< class T>
inline bool DerivativeObj<ThreeD,T>::IsPowerEqualTo(const int p1, const int p2, const int p3) const {
  return ( Power1==p1 && Power2==p2 && Power3==p3 );
}

// Friend functions
template< class T> inline
ostream& operator<<(ostream& os, const DerivativeObj<ThreeD,T>& Obj){

  os.width(4); os << Obj.P1(); 
  os.width(4); os << Obj.P2();
  os.width(4); os << Obj.P3();
  os.precision(15); os.width(25);

  return os << Obj.D();
}

template< class T> inline
void DerivativeObj<ThreeD,T>::Read_Derivative(istream &In_File){

  In_File.setf(ios::skipws);
  In_File >> Power1 >> Power2 >> Power3 >> ValueD; 
  In_File.unsetf(ios::skipws);
}

template< class T>
istream& operator>>(istream& os, DerivativeObj<ThreeD,T>& Obj){

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
const TaylorDerivativesContainer<ThreeD,T> operator+ (const TaylorDerivativesContainer<ThreeD,T>& left,
						    const TaylorDerivativesContainer<ThreeD,T>& right);

template<SpaceType ThreeD, class T>
const TaylorDerivativesContainer<ThreeD,T> operator- (const TaylorDerivativesContainer<ThreeD,T>& left,
						    const TaylorDerivativesContainer<ThreeD,T>& right);

template<SpaceType ThreeD, class T>
bool operator< (const TaylorDerivativesContainer<ThreeD,T>& left,
		const double& Value);

template<SpaceType ThreeD, class T>
bool operator<= (const TaylorDerivativesContainer<ThreeD,T>& left,
		 const double& Value);

template<SpaceType ThreeD, class T>
ostream & operator<< (ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj);

template<SpaceType ThreeD, class T>
istream & operator>> (istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj);

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
  short int ContainerSize;
  short int OrderOfRec;         /* OrderOfReconstruction */
  T phi;			/* Limiter value */
  T phi_copy;			/* Limiter value copied*/

  // Memory Management
  void allocate(const int NumberOfObjects);

public:
  // Default Constructor
  TaylorDerivativesContainer(void);
  // Overloaded constructor
  TaylorDerivativesContainer(const int OrderOfReconstruction);
  // Copy constructor
  TaylorDerivativesContainer(const TaylorDerivativesContainer<ThreeD,T> & rhs);
  // Assignment operator;
  TaylorDerivativesContainer<ThreeD,T> & operator=(const TaylorDerivativesContainer<ThreeD,T> & rhs);

  // GenerateContainer()
  void GenerateContainer(const int OrderOfReconstruction); /* Allocate memory for the derivatives and 
							      generate the power combinations */
  // free_memory()
  void free_memory(void); 		/* Delete the objects of the container */
  // Destructor
  ~TaylorDerivativesContainer(void){ free_memory();}

  // Characteristic parameters
  int size(void) const { return ContainerSize;}
  int FirstElem(void) const { return 0; }
  int LastElem(void) const { return ContainerSize-1;}
  int RecOrder(void) const { return OrderOfRec;}

  // IndexOrder(int,int,int) -> determines the position of the pointer to the element having the powers (p1,p2,p3)
  int IndexOrder(const int p1, const int p2, const int p3);
  const int IndexOrder(const int p1, const int p2, const int p3) const;

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ);
 
  // ComputeXGradientFor( )
  T ComputeXGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ); // RR: Not yet implemented
  // ComputeYGradientFor( )
  T ComputeYGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ); // RR: Not yet implemented
  // ComputeZGradientFor( )
  T ComputeZGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ); // RR: Not yet implemented

  // Reset limiter --> set limiter to ONE for all parameters
  void ResetLimiter(void){ phi.One();}
  void ResetFrozenLimiter(void){ phi_copy.One(); }
  
  /* Field access */
  const T & Limiter(void) const {return phi;}
  T & Limiter(void) {return phi;}
  const double & Limiter(const int Variable) const {return phi[Variable];}
  double & Limiter (const int Variable){return phi[Variable];}
  const T & Frozen_Limiter(void) const {return phi_copy;}
  T & Frozen_Limiter(void) {return phi_copy;}
  const double & Frozen_Limiter(const int Variable) const {return phi_copy[Variable];}
  double & Frozen_Limiter(const int Variable){return phi_copy[Variable];}
  void Make_Limiter_Copy(void){ phi_copy = phi; }
  void Make_Limiter_Copy(const int &Variable){ phi_copy[Variable] = phi[Variable]; }
  
  /* Overloaded Operators */
  Derivative & operator()(const int & position) {return DContainer[position];}
  const Derivative & operator()(const int & position) const {return DContainer[position];}

  T & operator()(const int p1, const int p2, const int p3);
  const T & operator()(const int p1, const int p2, const int p3) const;

  double & Value(const int position, const int parameter) { return DContainer[position].D(parameter); }
  const double & Value(const int position, const int parameter) const { return DContainer[position].D(parameter) };

  /* Friend functions */  
  friend const TaylorDerivativesContainer<ThreeD,T> operator+ <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
									  const TaylorDerivativesContainer<ThreeD,T>& right);
  
  friend const TaylorDerivativesContainer<ThreeD,T> operator- <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
									  const TaylorDerivativesContainer<ThreeD,T>& right);
  friend bool operator== <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
				     const TaylorDerivativesContainer<ThreeD,T>& right);
  
  friend bool operator!= <ThreeD,T> (const TaylorDerivativesContainer<ThreeD,T>& left,
				     const TaylorDerivativesContainer<ThreeD,T>& right);
  
  friend ostream & operator<< <ThreeD,T> (ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj);
  
  friend istream & operator>> <ThreeD,T> (istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj);
  
};

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer          *
 * Implementation of the Member Functions              *
 ******************************************************/

/* Default Constructor */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T>::TaylorDerivativesContainer(void):
  DContainer(NULL), ContainerSize(0), OrderOfRec(-1), phi(1.0), phi_copy(1.0){ }

/* Overloaded constructor */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T>::TaylorDerivativesContainer(const int OrderOfReconstruction):
  DContainer(NULL), ContainerSize(0), OrderOfRec(OrderOfReconstruction), phi(1.0), phi_copy(1.0)
{
  GenerateContainer(OrderOfRec);
}
 /* Copy constructor  */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<ThreeD,T> & rhs)
{

  // allocate memory for the new container
  allocate(rhs.size());
  
  // set the OrderOfRec
  OrderOfRec = rhs.RecOrder();

  /* Copy phi */
  phi = rhs.Limiter();
  phi_copy = rhs.phi_copy;

  // copy the values from the RHS
  for (int i=0; i<=rhs.LastElem(); ++i){
    DContainer[i] = rhs.DContainer[i];
  }
}

/* Allocate memory for the derivatives */
template<class T> inline
void TaylorDerivativesContainer<ThreeD,T>::allocate(const int NumberOfObjects)
{

  assert( NumberOfObjects >= 1 );
  if (ContainerSize != NumberOfObjects){
    /* free the memory if there is memory allocated */
    free_memory();
    
    ContainerSize = NumberOfObjects;
    /* create memory */    
    DContainer = new Derivative[ContainerSize];
  }

//  ContainerSize = NumberOfObjects;
//
//  if (ContainerSize > 0){
//    /* create memory */    
//    DContainer = new Derivative[ContainerSize];
//    
//    /* initialize state to ZERO */
//    Derivative ZERO_Derivative(0.0);
//    for(int i=0; i<=LastElem(); ++i){
//      DContainer[i] = ZERO_Derivative;
//    }
//  } else {
//    DContainer = NULL;
//  }

}

/* GenerateContainer --> initialize the derivatives and the power coefficients */
template<class T> inline
void TaylorDerivativesContainer<ThreeD,T>::GenerateContainer(const int OrderOfReconstruction)
{

  int NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)*(OrderOfReconstruction+3)/6;
  int p1,p2,p3,Position;
  T ZeroState(0.0);

  /* allocate memory for the derivatives */
  allocate(NumberOfDerivatives);

  // set the OrderOfRec
  OrderOfRec = (short int)OrderOfReconstruction;

  /* Set powers */
  for (p1=0,Position=0; p1<=OrderOfRec; ++p1)
    for (p2=0; p2<=OrderOfRec-p1; ++p2)
      for (p3=0; p3<=OrderOfRec-p1-p2; ++p3, ++Position){
	DContainer[Position].SetPowers(p1,p2,p3,true);
	DContainer[Position].SetValue(ZeroState);
    }
}

/* free_memory() */
template<class T> inline
void TaylorDerivativesContainer<ThreeD,T>::free_memory(void) {

  if (DContainer != NULL){ // Check if DContainer is not empty
    delete [] DContainer; DContainer = NULL;
    ContainerSize = 0;
    OrderOfRec = -1;
    phi = T(1.0);
    phi_copy = T(1.0);
  }
}

/* Assignment operator = */
template<class T> inline
TaylorDerivativesContainer<ThreeD,T> & 
TaylorDerivativesContainer<ThreeD,T>::operator=(const TaylorDerivativesContainer<ThreeD,T> & rhs){
  
  // !!! If the LHS container already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & rhs) return *this;

  // allocate memory if there isn't enough
  allocate(rhs.size());

  // set the OrderOfRec
  OrderOfRec = rhs.RecOrder();

  /* Copy phi */
  phi = rhs.Limiter();
  phi_copy = rhs.phi_copy;

  // copy the value from the RHS
  for (int i=0; i<=LastElem() ; ++i){
    DContainer[i] = rhs.DContainer[i];
  }

  return *this;
}

/* Access data */
/*IndexOrder(int,int,int) -> determines the position of the element having the power combination (p1,p2,p3) */
template<class T> inline
int TaylorDerivativesContainer<ThreeD,T>::IndexOrder(const int p1, const int p2, const int p3)
{
  int shift_p1 = 0;
  int shift_p2 = 0;
  int tmp_p2 = 0;
  int n, i;

  // Triangular numbers may be defined by f(n) = n*(n+1)/2 (or by the series f = 1,3,6,10.15,21,...)
  // This series of numbers is used to help shift the index of the derivatives container based on the p1 power.
  // Example: For OrderOfRec = 4 The shift_p1 value for p1=0 --> shift_p1 = 0         
  // -------                                        for p1=1 --> shift_p1 = 15        
  //                                                for p1=2 --> shift_p1 = 15+10     = 25
  //                                                for p1=3 --> shift_p1 = 15+10+6   = 31
  //                                                for p1=4 --> shift_p1 = 15+10+6+3 = 34

  if (p1+p2+p3 > OrderOfRec) {
    error_flag = 1; 
    cout << "\n " << CFFC_Version() 
	 << "Error in TaylorDerivativesContainer::IndexOrder -> Powers are out of bounds! ie. p1+p2+p3 > OrderOfRec, "
	 << "flag = " << error_flag << ".\n";
    return(error_flag);
  } /* endif */

  shift_p1 = 0;
  for (i=1; i<=p1; ++i){
    n = OrderOfRec-i+2;
    shift_p1 += n*(n+1)/2;
  }

  tmp_p2 = p1*p2;

  for ( i=2; i<=p2; ++i){
    tmp_p2 += (i-1);
  }
  shift_p2 = (OrderOfRec+1)*p2 - tmp_p2;

  // Return the position of the element
  return shift_p1 + shift_p2 + p3;
}


// operator(int) -> returns the derivative for which Power1 = p1 & Power2 = p2 & Power3 = p3
template<class T> inline
T & TaylorDerivativesContainer<ThreeD,T>::operator()(const int p1, const int p2, const int p3){

  return DContainer[IndexOrder(p1,p2,p3)].D();
}

// operator(int) -> returns the derivative for which p1_ = p1, p2_ = p2, p3_ = p3
template<class T> inline 
const T & TaylorDerivativesContainer<ThreeD,T>::operator()(const int p1, const int p2, const int p3) const{

  return DContainer[IndexOrder(p1,p2,p3)].D();
}


// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<class T> inline
T TaylorDerivativesContainer<ThreeD,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ){

  // initialize the solution state
  T Solution(0.0);

  int p1(0),p2(0),p3(0),Position(0);

  double DeltaXtoPower(1.0), DeltaYtoPower, DeltaZtoPower;

 for (p1=0,Position=0; p1<=OrderOfRec; ++p1){
    /* Reinitialize DeltaYtoPower */
    DeltaYtoPower = 1.0;

    for (p2=0; p2<=OrderOfRec-p1; ++p2){
      /* Reinitialize DeltaZtoPower */
      DeltaZtoPower = 1.0;
      
      for (p3=0; p3<=OrderOfRec-p1-p2; ++p3, ++Position){

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

 return (phi^Solution) + ((T(1.0)-phi)^DContainer[0].D());
}

// ComputeXGradientFor( ) :-> Compute the gradient in X direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<class T> inline
T TaylorDerivativesContainer<ThreeD,T>::ComputeXGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ){
  // RR: Not yet implemented
  return T(0.0);
}

// ComputeYGradientFor( ) :-> Compute the gradient in Y direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<class T> inline
T TaylorDerivativesContainer<ThreeD,T>::ComputeYGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ){
  // RR: Not yet implemented
  return T(0.0);
}

// ComputeZGradientFor( ) :-> Compute the gradient in ZX direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<class T> inline
T TaylorDerivativesContainer<ThreeD,T>::ComputeZGradientFor(const double DeltaX, const double DeltaY, const double DeltaZ){
  // RR: Not yet implemented
  return T(0.0);
}

// Friend functions
template<class T> inline
ostream & operator<< (ostream & os, const TaylorDerivativesContainer<ThreeD,T>& Obj){

  os.setf(ios::skipws);
  os.width(4);
  os << Obj.RecOrder() << endl;
  for(int i=0; i<=Obj.LastElem(); ++i)
    os << Obj(i) << endl;
  os.unsetf(ios::skipws);
  return os;
}

template<class T> inline
istream & operator>> (istream & os, TaylorDerivativesContainer<ThreeD,T>& Obj){

  int ReconstructionOrder;
  os.setf(ios::skipws);
  os >> ReconstructionOrder;
  /* Adjust the container size */
  Obj.GenerateContainer(ReconstructionOrder);
  for(int i=0; i<=Obj.LastElem(); ++i){
    os >> Obj(i);
  }
  os.unsetf(ios::skipws);
  return os;
}

// operator + -> summation of the derivatives
template<class T> inline
const TaylorDerivativesContainer<ThreeD,T> operator+(const TaylorDerivativesContainer<ThreeD,T>& left,
						   const TaylorDerivativesContainer<ThreeD,T>& right){

  TaylorDerivativesContainer<ThreeD,T> Temp(left.RecOrder());
  
  for (int i=0; i<=right.LastElem(); ++i){
    Temp(i).D() = left(i).D() + right(i).D();
  }
  
  return Temp;
}

// operator - -> difference of the derivatives
template<class T> inline
const TaylorDerivativesContainer<ThreeD,T> operator-(const TaylorDerivativesContainer<ThreeD,T>& left,
						   const TaylorDerivativesContainer<ThreeD,T>& right){
  
  TaylorDerivativesContainer<ThreeD,T> Temp(left.RecOrder());
  
  for (int i=0; i<=right.LastElem(); ++i){
    Temp(i).D() = left(i).D() - right(i).D();
  }
  
  return Temp;
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


// Specialization for "class T == double"

template<> inline
const double & DerivativeObj<ThreeD,double>::D(const int ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<ThreeD,double>::D(const int ) {
  return ValueD;
}

template<> inline
const double & TaylorDerivativesContainer<ThreeD,double>::Limiter(const int ) const {
  return phi;
}

template<> inline
double & TaylorDerivativesContainer<ThreeD,double>::Limiter(const int ) {
  return phi;
}

template<> inline
const double & TaylorDerivativesContainer<ThreeD,double>::Frozen_Limiter(const int ) const {
  return phi_copy;
}

template<> inline
double & TaylorDerivativesContainer<ThreeD,double>::Frozen_Limiter(const int ) {
  return phi_copy;
}

template<> inline
void TaylorDerivativesContainer<ThreeD,double>::ResetLimiter(void){
  phi = 1.0;
}

template<> inline
void TaylorDerivativesContainer<ThreeD,double>::ResetFrozenLimiter(void){
  phi_copy = 1.0;
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY,DeltaZ)
template<> inline
T TaylorDerivativesContainer<ThreeD,double>::ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ){

  // initialize the solution state
  double Solution(0.0);

  int p1(0),p2(0),p3(0),Position(0);

  double DeltaXtoPower(1.0), DeltaYtoPower, DeltaZtoPower;

 for (p1=0,Position=0; p1<=OrderOfRec; ++p1){
    /* Reinitialize DeltaYtoPower */
    DeltaYtoPower = 1.0;

    for (p2=0; p2<=OrderOfRec-p1; ++p2){
      /* Reinitialize DeltaZtoPower */
      DeltaZtoPower = 1.0;
      
      for (p3=0; p3<=OrderOfRec-p1-p2; ++p3, ++Position){

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

 return (phi*Solution) + ((1.0-phi)*DContainer[0].D());
}

#endif

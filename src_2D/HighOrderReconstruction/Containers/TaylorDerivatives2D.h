/*!\file TaylorDerivatives2D.h
  \brief Header file defining the specialization of 'TaylorDerivatives' template in 2D. 
  Obs: The factorial coefficients are incorporated in the values of the derivatives. */

#ifndef _TAYLORDERIVATIVES_2D_INCLUDED
#define _TAYLORDERIVATIVES_2D_INCLUDED

/* Include required C++ libraries. */
#include <cassert>

/* Using std namespace functions */
// None

/*****************************************************************************
// Define the partial specialization of the object class DerivativeObj in 2D *
*****************************************************************************/
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
 ostream& operator<< (ostream& os, const DerivativeObj<TwoD,T>& Obj);

template<SpaceType TwoD, class T >
istream& operator>> (istream& os, DerivativeObj<TwoD,T>& Obj);

/*******************************************************
 * CLASS Template: DerivativeObj                    *
 ******************************************************/
template< class T>
class DerivativeObj<TwoD,T>{

 private:
 short int Power1;  		/* the power coefficients */
 short int Power2;
 T ValueD;         		/* the Value of the Derivative, T = solution class */

 public:
 // Constructors
 DerivativeObj(void): Power1(0), Power2(0), ValueD() { };
 DerivativeObj(const double Val): Power1(0), Power2(0), ValueD(Val) { };
 DerivativeObj(const vector<int> & PPP_, const T & ValueD_);
 DerivativeObj(const int p1, const int p2, const T ValueD_): Power1(p1),Power2(p2),ValueD(ValueD_) { };
 // Copy constructor
 DerivativeObj( const DerivativeObj & rhs): Power1(rhs.Power1), Power2(rhs.Power2), ValueD(rhs.ValueD){ };
 ~DerivativeObj(){ };
 
 // Member functions
 void SetPowers(const int p1, const int p2) {Power1 = p1; Power2 = p2; }	/* set p1 & p2 */
 /* set p1 & p2 , passing by reference */
 void SetPowers(const int &p1, const int &p2, const bool) {Power1 = p1; Power2 = p2;}
 void SetPowers(const vector<int> & PPP_); /* set the powers */
 void SetValue(const T ValueD_) { ValueD = ValueD_; } /* set ValueD */
 void SetValue(const T & ValueD_, const bool ) { ValueD = ValueD_;} /* set ValueD, passing by reference */
 bool IsPowerEqualTo(const int p1, const int p2) const ; /* if the argument provided is equal to p1 & p2 returns true */

 void Read_Derivative(istream &In_File);

 // Access functions
 T & D(void) { return ValueD; } 		      /* return ValueD */
 const T & D(void) const {return ValueD; }
 const double & D(const int VarPosition) const { return ValueD[VarPosition]; }
 double & D(const int VarPosition){ return ValueD[VarPosition]; }
 short int & P1(void) { return Power1; }	/* returns p1 */
 const short int & P1(void) const {return Power1;}
 short int & P2(void) { return Power2; }	/* returns p2 */
 const short int & P2(void) const {return Power2;}

 /* Overloaded operators */
 // Assignment operator
 DerivativeObj& operator=(const DerivativeObj<TwoD,T>& rhs);
 
 // Friend functions & operators
 friend bool operator== <TwoD,T> (const DerivativeObj<TwoD,T>& left,
				  const DerivativeObj<TwoD,T>& right);
 friend bool operator!= <TwoD,T> (const DerivativeObj<TwoD,T>& left,
				  const DerivativeObj<TwoD,T>& right);
 friend ostream& operator<< <TwoD,T> (ostream& os, const DerivativeObj<TwoD,T>& Obj);

 friend istream& operator>> <TwoD,T> (istream& os, DerivativeObj<TwoD,T>& Obj);
};

// CLASS DerivativeObj

// Constructor
template <class T> inline
DerivativeObj<TwoD,T>::DerivativeObj(const vector<int> & PPP_, const T & ValueD_){
   require(PPP_.size() == 2,"Constructor DerivativeObj failed due to difference in dimensions.\n");
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
inline void DerivativeObj<TwoD,T>::SetPowers(const vector<int> & PPP_){
  require(PPP_.size() >= 1, "DerivativeObj::SetPowers(vector<int>) failed. Not enough elements.\n");
  Power1 = PPP_[0];
  Power2 = PPP_[1];
}

// IsPowerEqualTo()
template< class T>
inline bool DerivativeObj<TwoD,T>::IsPowerEqualTo(const int p1, const int p2) const {
  return ( Power1==p1 && Power2==p2 );
}

// Friend functions
template< class T> inline
ostream& operator<<(ostream& os, const DerivativeObj<TwoD,T>& Obj){

  os.width(4); os << Obj.P1(); 
  os.width(4); os << Obj.P2();
  os.precision(15); os.width(25);

  return os << Obj.D();
}

template< class T> inline
void DerivativeObj<TwoD,T>::Read_Derivative(istream &In_File){

  In_File.setf(ios::skipws);
  In_File >> Power1 >> Power2 >> ValueD; 
  In_File.unsetf(ios::skipws);
}

template< class T>
istream& operator>>(istream& os, DerivativeObj<TwoD,T>& Obj){

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
const TaylorDerivativesContainer<TwoD,T> operator+ (const TaylorDerivativesContainer<TwoD,T>& left,
						    const TaylorDerivativesContainer<TwoD,T>& right);

template<SpaceType TwoD, class T>
const TaylorDerivativesContainer<TwoD,T> operator- (const TaylorDerivativesContainer<TwoD,T>& left,
						    const TaylorDerivativesContainer<TwoD,T>& right);

template<SpaceType TwoD, class T>
bool operator< (const TaylorDerivativesContainer<TwoD,T>& left,
		const double& Value);

template<SpaceType TwoD, class T>
bool operator<= (const TaylorDerivativesContainer<TwoD,T>& left,
		 const double& Value);

template<SpaceType TwoD, class T>
ostream & operator<< (ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj);

template<SpaceType TwoD, class T>
istream & operator>> (istream & os, TaylorDerivativesContainer<TwoD,T>& Obj);

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
  short int ContainerSize;
  short int OrderOfRec;  	/* OrderOfReconstruction */
  T phi;			/* Limiter value */
  T phi_copy;			/* Limiter value copied */

  // Memory Management
  void allocate(const int NumberOfObjects);

 public:
  // Default Constructor
  TaylorDerivativesContainer(void);
  // Overloaded constructor
  TaylorDerivativesContainer(const int OrderOfReconstruction);
  // Copy constructor
  TaylorDerivativesContainer(const TaylorDerivativesContainer<TwoD,T> & rhs);
  // Assignment operator;
  TaylorDerivativesContainer<TwoD,T> & operator=(const TaylorDerivativesContainer<TwoD,T> & rhs);

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

  // IndexOrder(int,int) -> determines the position of the pointer to the element having the powers (p1,p2)
  int IndexOrder(const int p1, const int p2) const;

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX, const double DeltaY);

  // ComputeXGradientFor( )
  T ComputeXGradientFor(const double DeltaX, const double DeltaY);

  // ComputeYGradientFor( )
  T ComputeYGradientFor(const double DeltaX, const double DeltaY);

  //! Compute X-gradient for a specific solution variable
  double ComputeXGradientFor(const double & DeltaX, const double & DeltaY, const unsigned & parameter) const;

  //! Compute Y-gradient for a specific solution variable
  double ComputeYGradientFor(const double & DeltaX, const double & DeltaY, const unsigned & parameter) const;

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
  Derivative & operator()(const int position) {return DContainer[position];}
  const Derivative & operator()(const int position) const {return DContainer[position];}

  T & operator()(const int p1, const int p2);
  const T & operator()(const int p1, const int p2) const;

  double & Value(const int position, const int parameter) { return DContainer[position].D(parameter); }
  const double & Value(const int position, const int parameter) const { return DContainer[position].D(parameter); }

  /* Friend functions */
  friend const TaylorDerivativesContainer<TwoD,T> operator+ <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
								      const TaylorDerivativesContainer<TwoD,T>& right);

  friend const TaylorDerivativesContainer<TwoD,T> operator- <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
								      const TaylorDerivativesContainer<TwoD,T>& right);

  friend bool operator== <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
				   const TaylorDerivativesContainer<TwoD,T>& right);

  friend bool operator!= <TwoD,T> (const TaylorDerivativesContainer<TwoD,T>& left,
				   const TaylorDerivativesContainer<TwoD,T>& right);

  friend ostream & operator<< <TwoD,T> (ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj);

  friend istream & operator>> <TwoD,T> (istream & os, TaylorDerivativesContainer<TwoD,T>& Obj);

};

/*******************************************************
 * CLASS Template: TaylorDerivativesContainer          *
 * Implementation of the Member Functions              *
 ******************************************************/
/* Default Constructor */
template<class T> inline
TaylorDerivativesContainer<TwoD,T>::TaylorDerivativesContainer(void):
  DContainer(NULL), ContainerSize(0), OrderOfRec(-1), phi(1.0), phi_copy(1.0){ }

/* Overloaded constructor */
template<class T> inline
TaylorDerivativesContainer<TwoD,T>::TaylorDerivativesContainer(const int OrderOfReconstruction):
  DContainer(NULL), ContainerSize(0), OrderOfRec(OrderOfReconstruction), phi(1.0), phi_copy(1.0)
{
  GenerateContainer(OrderOfRec);
}

/* Copy constructor  */
template<class T> inline
TaylorDerivativesContainer<TwoD,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<TwoD,T> & rhs)
  :DContainer(NULL), ContainerSize(0), OrderOfRec(-1), phi(1.0), phi_copy(1.0) {

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
void TaylorDerivativesContainer<TwoD,T>::allocate(const int NumberOfObjects)
{
  assert( NumberOfObjects >= 1 );
  if (ContainerSize != NumberOfObjects){
    /* free the memory if there is memory allocated */
    free_memory();
    
    ContainerSize = NumberOfObjects;
    /* create memory */    
    DContainer = new Derivative[ContainerSize];
  }
}

/* GenerateContainer --> initialize the derivatives and the power coefficients */
template<class T> inline
void TaylorDerivativesContainer<TwoD,T>::GenerateContainer(const int OrderOfReconstruction)
{

  int NumberOfDerivatives((OrderOfReconstruction+1)*(OrderOfReconstruction+2)/2);
  int p1,p2,Position;
  T ZeroState(0.0);

  /* allocate memory for the derivatives */
  allocate(NumberOfDerivatives); 

  // set the OrderOfRec
  OrderOfRec = (short int)OrderOfReconstruction;

  /* Set powers */
  for (p1=0,Position=0; p1<=OrderOfRec; ++p1){
    for (p2=0; p2<=OrderOfRec-p1; ++p2, ++Position){
      DContainer[Position].SetPowers(p1,p2,true);
      DContainer[Position].SetValue(ZeroState);
    }
  }
}

/* free_memory() */
template<class T> inline
void TaylorDerivativesContainer<TwoD,T>::free_memory(void) {

  if (DContainer != NULL){ // Check if DContainer is empty
    delete [] DContainer; DContainer = NULL;
    ContainerSize = 0;
    OrderOfRec = -1;
    phi = T(1.0);
    phi_copy = T(1.0);
  }
}

/* Assignment operator = */
template<class T> inline
TaylorDerivativesContainer<TwoD,T> & 
TaylorDerivativesContainer<TwoD,T>::operator=(const TaylorDerivativesContainer<TwoD,T> & rhs){

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
/*IndexOrder(int,int) -> determines the position of the element having the power combination (p1,p2) */
template<class T> inline
  int TaylorDerivativesContainer<TwoD,T>::IndexOrder(const int p1, const int p2) const
{
  int sum(0), n;

  for (n=2; n<=p1; ++n){
    sum += (n - 1);
  }
  return (OrderOfRec + 1)*p1  - sum + p2;
}

// operator(int) -> returns the derivative for which Power1 = p1 & Power2 = p2
template<class T> inline
T & TaylorDerivativesContainer<TwoD,T>::operator()(const int p1, const int p2){

  return DContainer[IndexOrder(p1,p2)].D();
}

// operator(int) -> returns the derivative for which Power1 = p1 & Power2 = p2
template<class T> inline 
const T & TaylorDerivativesContainer<TwoD,T>::operator()(const int p1, const int p2) const{

  return DContainer[IndexOrder(p1,p2)].D();
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY)
template<class T> inline
T TaylorDerivativesContainer<TwoD,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY){

  double DeltaXSquare, DeltaYSquare;


  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;
    
    return ( DContainer[0].D() + (phi^(DeltaY*DContainer[1].D() + DeltaYSquare*DContainer[2].D()
				       + DeltaYSquare*DeltaY*DContainer[3].D() + DeltaYSquare*DeltaYSquare*DContainer[4].D()
				       + DeltaX*DContainer[5].D() + DeltaX*DeltaY*DContainer[6].D()
				       + DeltaX*DeltaYSquare*DContainer[7].D() + DeltaX*DeltaYSquare*DeltaY*DContainer[8].D()
				       + DeltaXSquare*DContainer[9].D() + DeltaXSquare*DeltaY*DContainer[10].D()
				       + DeltaXSquare*DeltaYSquare*DContainer[11].D() + DeltaX*DeltaXSquare*DContainer[12].D() 
				       + DeltaX*DeltaXSquare*DeltaY*DContainer[13].D()
				       + DeltaXSquare*DeltaXSquare*DContainer[14].D() ) ) );
    break;
    

  case 3: 			/* Third Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[0].D() + (phi^(DeltaY*DContainer[1].D() + DeltaYSquare*DContainer[2].D()
				       + DeltaYSquare*DeltaY*DContainer[3].D() + DeltaX*DContainer[4].D()
				       + DeltaX*DeltaY*DContainer[5].D()
				       + DeltaX*DeltaYSquare*DContainer[6].D() + DeltaXSquare*DContainer[7].D() 
				       + DeltaXSquare*DeltaY*DContainer[8].D()
				       + DeltaX*DeltaXSquare*DContainer[9].D()) ) );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[0].D() + (phi^ (DeltaY*DContainer[1].D() + DeltaY*DeltaY*DContainer[2].D() +
					DeltaX*DContainer[3].D() + DeltaX*DeltaY*DContainer[4].D() +
					DeltaX*DeltaX*DContainer[5].D()) ) );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[0].D() + (phi^(DeltaY*DContainer[1].D() + DeltaX*DContainer[2].D())) );
    break;

  case 0:
    return DContainer[0].D();
    break;

  default:
    /* Orders higher than FOUR */

    // initialize the solution state
    T Solution(0.0);
    double DeltaXtoPower(1.0), DeltaYtoPower;
    int p1,p2,Position;
    for (p1=0,Position=0; p1<=OrderOfRec; ++p1){
      /* Reinitialize DeltaYtoPower */
      DeltaYtoPower = 1.0;
      for (p2=0; p2<=OrderOfRec-p1; ++p2, ++Position){
	/* Update solution */
	Solution += DeltaXtoPower*DeltaYtoPower*DContainer[Position].D();

	/* Update DeltaYtoPower */
	DeltaYtoPower *= DeltaY;
      }
      /* Update DeltaXtoPower */
      DeltaXtoPower *= DeltaX;
    }

    return (phi^Solution) + ((T(1.0)-phi)^DContainer[0].D());

  }
}

// ComputeXGradientFor( ) :-> Compute the gradient in X direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY)
template<class T> inline
T TaylorDerivativesContainer<TwoD,T>::ComputeXGradientFor(const double DeltaX, const double DeltaY){

  double DeltaXSquare, DeltaYSquare;

  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[5].D() + DeltaY*DContainer[6].D()
	     + DeltaYSquare*DContainer[7].D() + DeltaY*DeltaYSquare*DContainer[8].D() 
	     + 2.0*DeltaX*DContainer[9].D() + 2.0*DeltaX*DeltaY*DContainer[10].D() 
	     + 2.0*DeltaX*DeltaYSquare*DContainer[11].D() + 3.0*DeltaXSquare*DContainer[12].D()
	     + 3.0*DeltaXSquare*DeltaY*DContainer[13].D() + 4.0*DeltaXSquare*DeltaX*DContainer[14].D() );
    break;

  case 3: 			/* Third Order Reconstruction */
    return ( DContainer[4].D() + DeltaY*DContainer[5].D()
	     + DeltaY*DeltaY*DContainer[6].D() + 2.0*DeltaX*DContainer[7].D() 
	     + 2.0*DeltaX*DeltaY*DContainer[8].D() + 3.0*DeltaX*DeltaX*DContainer[9].D() );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[3].D() + DeltaY*DContainer[4].D() + 2.0*DeltaX*DContainer[5].D() );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[2].D() );
    break;

  case 0:
    return T(0.0);
    break;

  default:
    /* Orders higher than FOUR */
    return T(0.0);			/* TO DO LATER! */

  }
}

// ComputeYGradientFor( ) :-> Compute the gradient in Y direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY)
template<class T> inline
T TaylorDerivativesContainer<TwoD,T>::ComputeYGradientFor(const double DeltaX, const double DeltaY){

  double DeltaXSquare, DeltaYSquare;

  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[1].D() + 2.0*DeltaY*DContainer[2].D()
	     + 3.0*DeltaYSquare*DContainer[3].D() + 4.0*DeltaYSquare*DeltaY*DContainer[4].D()
	     + DeltaX*DContainer[6].D() + 2.0*DeltaX*DeltaY*DContainer[7].D() 
	     + 3.0*DeltaX*DeltaYSquare*DContainer[8].D() + DeltaXSquare*DContainer[10].D() 
	     + 2.0*DeltaXSquare*DeltaY*DContainer[11].D() + DeltaX*DeltaXSquare*DContainer[13].D() );
    break;

  case 3: 			/* Third Order Reconstruction */
    return ( DContainer[1].D() + 2.0*DeltaY*DContainer[2].D()
	     + 3.0*DeltaY*DeltaY*DContainer[3].D() + DeltaX*DContainer[5].D()
	     + 2.0*DeltaX*DeltaY*DContainer[6].D() + DeltaX*DeltaX*DContainer[8].D() );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[1].D() + 2.0*DeltaY*DContainer[2].D() + DeltaX*DContainer[4].D() );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[1].D() );
    break;

  case 0:
    return T(0.0);
    break;

  default:
    /* Orders higher than FOUR */
    return T(0.0);			/* TO DO LATER! */
  }
}

// ComputeXGradientFor( ) :-> Compute the gradient in X direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY) and a specific solution variable
template<class T> inline
double TaylorDerivativesContainer<TwoD,T>::ComputeXGradientFor(const double &DeltaX, const double &DeltaY,
							       const unsigned & parameter) const {

  double DeltaXSquare, DeltaYSquare;

  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[5].D(parameter) + DeltaY*DContainer[6].D(parameter)
	     + DeltaYSquare*DContainer[7].D(parameter) + DeltaY*DeltaYSquare*DContainer[8].D(parameter) 
	     + 2.0*DeltaX*DContainer[9].D(parameter) + 2.0*DeltaX*DeltaY*DContainer[10].D(parameter) 
	     + 2.0*DeltaX*DeltaYSquare*DContainer[11].D(parameter) + 3.0*DeltaXSquare*DContainer[12].D(parameter)
	     + 3.0*DeltaXSquare*DeltaY*DContainer[13].D(parameter) + 4.0*DeltaXSquare*DeltaX*DContainer[14].D(parameter) );
    break;

  case 3: 			/* Third Order Reconstruction */
    return ( DContainer[4].D(parameter) + DeltaY*DContainer[5].D(parameter)
	     + DeltaY*DeltaY*DContainer[6].D(parameter) + 2.0*DeltaX*DContainer[7].D(parameter) 
	     + 2.0*DeltaX*DeltaY*DContainer[8].D(parameter) + 3.0*DeltaX*DeltaX*DContainer[9].D(parameter) );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[3].D(parameter) + DeltaY*DContainer[4].D(parameter) + 2.0*DeltaX*DContainer[5].D(parameter) );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[2].D(parameter) );
    break;

  case 0:
    return 0.0;
    break;

  default:
    /* Orders higher than FOUR */
    return 0.0;			/* TO DO LATER! */

  }
}

// ComputeYGradientFor( ) :-> Compute the gradient in Y direction of
// the Taylor series expansion for a particular distance (DeltaX,DeltaY) and a specific solution variable
template<class T> inline
double TaylorDerivativesContainer<TwoD,T>::ComputeYGradientFor(const double &DeltaX, const double &DeltaY,
							       const unsigned & parameter) const {

  double DeltaXSquare, DeltaYSquare;

  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[1].D(parameter) + 2.0*DeltaY*DContainer[2].D(parameter)
	     + 3.0*DeltaYSquare*DContainer[3].D(parameter) + 4.0*DeltaYSquare*DeltaY*DContainer[4].D(parameter)
	     + DeltaX*DContainer[6].D(parameter) + 2.0*DeltaX*DeltaY*DContainer[7].D(parameter) 
	     + 3.0*DeltaX*DeltaYSquare*DContainer[8].D(parameter) + DeltaXSquare*DContainer[10].D(parameter) 
	     + 2.0*DeltaXSquare*DeltaY*DContainer[11].D(parameter) + DeltaX*DeltaXSquare*DContainer[13].D(parameter) );
    break;

  case 3: 			/* Third Order Reconstruction */
    return ( DContainer[1].D(parameter) + 2.0*DeltaY*DContainer[2].D(parameter)
	     + 3.0*DeltaY*DeltaY*DContainer[3].D(parameter) + DeltaX*DContainer[5].D(parameter)
	     + 2.0*DeltaX*DeltaY*DContainer[6].D(parameter) + DeltaX*DeltaX*DContainer[8].D(parameter) );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[1].D(parameter) + 2.0*DeltaY*DContainer[2].D(parameter) + DeltaX*DContainer[4].D(parameter) );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[1].D(parameter) );
    break;

  case 0:
    return 0.0;
    break;

  default:
    /* Orders higher than FOUR */
    return 0.0;			/* TO DO LATER! */
  }
}


// Friend functions
template<class T> inline
ostream & operator<< (ostream & os, const TaylorDerivativesContainer<TwoD,T>& Obj){

  os.setf(ios::skipws);
  os.width(4);
  os << Obj.RecOrder() << endl;
  for(int i=0; i<=Obj.LastElem(); ++i)
    os << Obj(i) << endl;
  os.unsetf(ios::skipws);
  return os;
}

template<class T> inline
istream & operator>> (istream & os, TaylorDerivativesContainer<TwoD,T>& Obj){

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
const TaylorDerivativesContainer<TwoD,T> operator+(const TaylorDerivativesContainer<TwoD,T>& left,
						   const TaylorDerivativesContainer<TwoD,T>& right){

  TaylorDerivativesContainer<TwoD,T> Temp(left.RecOrder());
  
  for (int i=0; i<=right.LastElem(); ++i){
    Temp(i).D() = left(i).D() + right(i).D();
  }
  
  return Temp;
}

// operator - -> difference of the derivatives
template<class T> inline
const TaylorDerivativesContainer<TwoD,T> operator-(const TaylorDerivativesContainer<TwoD,T>& left,
						   const TaylorDerivativesContainer<TwoD,T>& right){
  
  TaylorDerivativesContainer<TwoD,T> Temp(left.RecOrder());
  
  for (int i=0; i<=right.LastElem(); ++i){
    Temp(i).D() = left(i).D() - right(i).D();
  }
  
  return Temp;
}

// operator ==
template<class T> inline
bool operator==(const TaylorDerivativesContainer<TwoD,T>& left,
		const TaylorDerivativesContainer<TwoD,T>& right){

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
inline bool operator!=(const TaylorDerivativesContainer<TwoD,T>& left,
		       const TaylorDerivativesContainer<TwoD,T>& right){

  return !(left == right);
}


// Specialization for "class T == double"
template<> inline
const double & DerivativeObj<TwoD,double>::D(const int ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<TwoD,double>::D(const int ) {
  return ValueD;
}

template<> inline
const double & TaylorDerivativesContainer<TwoD,double>::Limiter(const int ) const {
  return phi;
}

template<> inline
double & TaylorDerivativesContainer<TwoD,double>::Limiter(const int ) {
  return phi;
}

template<> inline
const double & TaylorDerivativesContainer<TwoD,double>::Frozen_Limiter(const int ) const {
  return phi_copy;
}

template<> inline
double & TaylorDerivativesContainer<TwoD,double>::Frozen_Limiter(const int ) {
  return phi_copy;
}

template<> inline
void TaylorDerivativesContainer<TwoD,double>::ResetLimiter(void){
  phi = 1.0;
}

template<> inline
void TaylorDerivativesContainer<TwoD,double>::ResetFrozenLimiter(void){
  phi_copy = 1.0;
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance (DeltaX,DeltaY)
template<> inline
double TaylorDerivativesContainer<TwoD,double>::ComputeSolutionFor(const double DeltaX, const double DeltaY){

  double DeltaXSquare, DeltaYSquare;


  switch(OrderOfRec){
  case 4:			/* Fourth Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[0].D() + (phi*(DeltaY*DContainer[1].D() + DeltaYSquare*DContainer[2].D()
				       + DeltaYSquare*DeltaY*DContainer[3].D() + DeltaYSquare*DeltaYSquare*DContainer[4].D()
				       + DeltaX*DContainer[5].D() + DeltaX*DeltaY*DContainer[6].D()
				       + DeltaX*DeltaYSquare*DContainer[7].D() + DeltaX*DeltaYSquare*DeltaY*DContainer[8].D()
				       + DeltaXSquare*DContainer[9].D() + DeltaXSquare*DeltaY*DContainer[10].D()
				       + DeltaXSquare*DeltaYSquare*DContainer[11].D() + DeltaX*DeltaXSquare*DContainer[12].D() 
				       + DeltaX*DeltaXSquare*DeltaY*DContainer[13].D()
				       + DeltaXSquare*DeltaXSquare*DContainer[14].D() ) ) );
    break;

  case 3: 			/* Third Order Reconstruction */
    DeltaXSquare = DeltaX*DeltaX;
    DeltaYSquare = DeltaY*DeltaY;

    return ( DContainer[0].D() + (phi*(DeltaY*DContainer[1].D() + DeltaYSquare*DContainer[2].D()
				       + DeltaYSquare*DeltaY*DContainer[3].D() + DeltaX*DContainer[4].D()
				       + DeltaX*DeltaY*DContainer[5].D()
				       + DeltaX*DeltaYSquare*DContainer[6].D() + DeltaXSquare*DContainer[7].D() 
				       + DeltaXSquare*DeltaY*DContainer[8].D()
				       + DeltaX*DeltaXSquare*DContainer[9].D()) ) );
    break;

  case 2:			/* Second Order Reconstruction */
    
    return ( DContainer[0].D() + (phi* (DeltaY*DContainer[1].D() + DeltaY*DeltaY*DContainer[2].D() +
					DeltaX*DContainer[3].D() + DeltaX*DeltaY*DContainer[4].D() +
					DeltaX*DeltaX*DContainer[5].D()) ) );
    break;

  case 1:			/* First Order Reconstruction */

    return ( DContainer[0].D() + (phi*(DeltaY*DContainer[1].D() + DeltaX*DContainer[2].D())) );
    break;

  case 0:
    return DContainer[0].D();
    break;

  default:
    /* Orders higher than THREE */

    // initialize the solution state
    double Solution(0.0);
    double DeltaXtoPower(1.0), DeltaYtoPower;
    int p1,p2,Position;
    for (p1=0,Position=0; p1<=OrderOfRec; ++p1){
      /* Reinitialize DeltaYtoPower */
      DeltaYtoPower = 1.0;
      for (p2=0; p2<=OrderOfRec-p1; ++p2, ++Position){
	/* Update solution */
	Solution += DeltaXtoPower*DeltaYtoPower*DContainer[Position].D();

	/* Update DeltaYtoPower */
	DeltaYtoPower *= DeltaY;
      }
      /* Update DeltaXtoPower */
      DeltaXtoPower *= DeltaX;
    }

    return (phi*Solution) + ((1.0 - phi)*DContainer[0].D());

  }
}

#endif

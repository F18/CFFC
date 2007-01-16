/* TaylorDerivatives.h: Header file defining the array of derivatives for a Taylor series expansion.
   Obs: The factorial coefficients are incorporated in the value of the derivative*/

#ifndef _TAYLORDERIVATIVES_TEMPLATE_INCLUDED
#define _TAYLORDERIVATIVES_TEMPLATE_INCLUDED

#include <iostream>
#include <vector>

#include "include/require.h"
#include "include/TypeDefinition.h"

// Define the object class: DerivativeObj
template< SpaceType SpaceDimension, class T>
class DerivativeObj;

/************************************************
*     Friend Functions : DerivativeObj          *
************************************************/
template< SpaceType SpaceDimension, class T >
bool operator==(const DerivativeObj<SpaceDimension,T>& left, const DerivativeObj<SpaceDimension,T>& right);

template< SpaceType SpaceDimension, class T >
bool operator!=(const DerivativeObj<SpaceDimension,T>& left, const DerivativeObj<SpaceDimension,T>& right);

template< SpaceType SpaceDimension, class T >
std::ostream& operator<< (std::ostream& os, const DerivativeObj<SpaceDimension,T>& Obj);

/*******************************************************
 * TEMPLETIZED CLASS: DerivativeObj                    *
 ******************************************************/
template< SpaceType SpaceDimension, class T>
class DerivativeObj{

 private:
 std::vector<int> PPP;		/* the power coefficients */
 T ValueD;         		/* the Value of the Derivative, T = solution class */

 public:
 // Constructors
 DerivativeObj(): ValueD() { allocate();}
 DerivativeObj(const double & Val): ValueD(Val) { allocate();}
 DerivativeObj(const std::vector<int> PPP_, const T ValueD_){
#ifndef __No_Checking__
   require(PPP_.size() == SpaceDimension,"Constructor DerivativeObj failed due to difference in dimensions.\n");
#endif
   PPP = PPP_;
   ValueD = ValueD_;
 };
 // Copy constructor
 DerivativeObj( const DerivativeObj & rhs): PPP(rhs.PPP), ValueD(rhs.ValueD){ };
 ~DerivativeObj(){ };
 
 // Member functions
 void allocate();		/* allocate space for the object */
 void SetPowers(const unsigned p1);	/* set p1 */
 void SetPowers(const unsigned p1, const unsigned p2); /* set p1 and p2 */
 void SetPowers(const unsigned p1, const unsigned p2, const unsigned p3); /* Set p1, p2 and p3*/
 void SetPowers(const std::vector<unsigned> PPP_); /* set the powers */
 void SetValue(const T ValueD_) { ValueD = ValueD_; } /* Set ValueD */
 T & D( ) { return ValueD; } 		/* return ValueD */
 const double & D(const unsigned VarPosition) const { return ValueD[VarPosition]; }
 double & D(const unsigned VarPosition){ return ValueD[VarPosition]; }

 unsigned P1(void) const ;		/* returns p1 */
 unsigned P2(void) const ;		/* returns p2 */
 unsigned P3(void) const ;		/* returns p3 */
 // std::vector<int> & Powers( ) const;
 bool IsPowerEqualTo(const unsigned p1) const ; /* if the argument provided is equal to p1 returns true */
 bool IsPowerEqualTo(const unsigned p1, const unsigned p2 ) const ;
 bool IsPowerEqualTo(const unsigned p1, const unsigned p2 , const unsigned p3) const ;


  /* Overloaded operators */
  // Assignment operator
 DerivativeObj& operator=(const DerivativeObj<SpaceDimension,T>& right);
 
 // equality operator
 friend bool operator== <SpaceDimension,T> (const DerivativeObj<SpaceDimension,T>& left,
					    const DerivativeObj<SpaceDimension,T>& right);
 friend bool operator!= <SpaceDimension,T> (const DerivativeObj<SpaceDimension,T>& left,
					    const DerivativeObj<SpaceDimension,T>& right);

 /*  ostream operator */
 friend std::ostream& operator<< <SpaceDimension,T> (std::ostream& os, const DerivativeObj<SpaceDimension,T>& Obj);
};

// CLASS DerivativeObj
// Assignment operator
template<SpaceType SpaceDimension, class T>
inline DerivativeObj<SpaceDimension,T>&
DerivativeObj<SpaceDimension,T>::operator=(const DerivativeObj<SpaceDimension,T>& rhs){
  if(this == &rhs) return *this;
  PPP = rhs.PPP;
  ValueD = rhs.ValueD;
  return *this;
}

// allocate()
template<SpaceType SpaceDimension, class T>
inline void DerivativeObj<SpaceDimension,T>::allocate(){
  switch(SpaceDimension){
  case OneD:
    PPP.assign(1,0);
    break;
  case TwoD:
    PPP.assign(2,0);
    break;
  case ThreeD:
    PPP.assign(3,0);
    break;
  }
}

// SetPowers()
template<SpaceType SpaceDimension, class T>
inline void DerivativeObj<SpaceDimension,T>::SetPowers(const std::vector<unsigned> PPP_){
#ifndef __No_Checking__
  require(SpaceDimension == PPP_.size(), "DerivativeObj::SetPowers(vector<int>) failed. The space dimension is different\n");
#endif
  PPP = PPP_;
}

template<SpaceType SpaceDimension, class T>
inline void DerivativeObj<SpaceDimension,T>::SetPowers(const unsigned p1){
#ifndef __No_Checking__
  require(SpaceDimension == OneD, "DerivativeObj::SetPowers(int) failed. The space dimension is different\n");
#endif
  PPP[0] = p1;
}

template<SpaceType SpaceDimension, class T>
inline void DerivativeObj<SpaceDimension,T>::SetPowers(const unsigned p1, const unsigned p2){
#ifndef __No_Checking__
  require(SpaceDimension == TwoD, "DerivativeObj::SetPowers(unsigned,unsigned) failed. The space dimension is different\n");
#endif
  PPP[0] = p1;
  PPP[1] = p2;
}

template<SpaceType SpaceDimension, class T>
inline void DerivativeObj<SpaceDimension,T>::SetPowers(const unsigned p1, const unsigned p2, const unsigned p3){
#ifndef __No_Checking__
  require(SpaceDimension == ThreeD, "DerivativeObj::SetPowers(unsigned,unsigned,unsigned) failed. The space dimension is different\n");
#endif
  PPP[0] = p1;
  PPP[1] = p2;
  PPP[2] = p3;
}

// IsPowerEqualTo
template<SpaceType SpaceDimension, class T>
inline bool DerivativeObj<SpaceDimension,T>::IsPowerEqualTo(const unsigned p1) const {
#ifndef __No_Checking__
  require(SpaceDimension == OneD, "DerivativeObj::IsPowerEqualTo(unsigned) failed. The space dimension is different\n");
#endif
  return PPP[0]==p1;
}

template<SpaceType SpaceDimension, class T>
inline bool DerivativeObj<SpaceDimension,T>::IsPowerEqualTo(const unsigned p1, const unsigned p2) const {
#ifndef __No_Checking__
  require(SpaceDimension == TwoD, "DerivativeObj::IsPowerEqualTo(unsigned,unsigned) failed. The space dimension is different\n");
#endif
  return (PPP[0]==p1)&&(PPP[1]==p2);
}

template<SpaceType SpaceDimension, class T>
inline bool DerivativeObj<SpaceDimension,T>::IsPowerEqualTo(const unsigned p1, const unsigned p2, const unsigned p3) const {
#ifndef __No_Checking__
  require(SpaceDimension == ThreeD, "DerivativeObj::IsPowerEqualTo(unsigned,unsigned,unsigned) failed. The space dimension is different\n");
#endif
  return (PPP[0]==p1)&&(PPP[1]==p2)&&(PPP[2]==p3);
}

// P1()
template<SpaceType SpaceDimension, class T>
inline unsigned DerivativeObj<SpaceDimension,T>::P1(void) const {
  return PPP[0];
}

// P2()
template<SpaceType SpaceDimension, class T>
inline unsigned DerivativeObj<SpaceDimension,T>::P2(void) const {
#ifndef __No_Checking__
  require(SpaceDimension != OneD,"DerivativeObj::P2( ) failed. The space dimension is different\n");
#endif
  return PPP[1];
}

// P3()
template<SpaceType SpaceDimension, class T>
inline unsigned DerivativeObj<SpaceDimension,T>::P3(void) const {
#ifndef __No_Checking__
  require(SpaceDimension == ThreeD,"DerivativeObj::P3( ) failed. The space dimension is different\n");
#endif
  return PPP[2];
}


// Friend functions
template<SpaceType SpaceDimension, class T> inline
std::ostream& operator<<(std::ostream& os, const DerivativeObj<SpaceDimension,T>& Obj){
  os << "powers= ";
  for (int i=0; i<Obj.PPP.size(); ++i)
    os << Obj.PPP[i] << " ,";
  return os << "\tValueD=" << Obj.ValueD << std::endl;
}

// operator ==
template<SpaceType SpaceDimension, class T> inline
bool operator==(const DerivativeObj<SpaceDimension,T>& left,
		const DerivativeObj<SpaceDimension,T>& right){
  return (left.PPP==right.PPP)&&(left.ValueD==right.ValueD);
}

// operator !=
template<SpaceType SpaceDimension, class T> inline
bool operator!=(const DerivativeObj<SpaceDimension,T>& left,
		const DerivativeObj<SpaceDimension,T>& right){
  return !(left == right);
}

template<> inline
const double & DerivativeObj<ThreeD,double>::D(const unsigned ) const {
  return ValueD;
}

template<> inline
double & DerivativeObj<ThreeD,double>::D(const unsigned ) {
  return ValueD;
}

/*******************************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer Declaration       *
 ******************************************************************/
template<SpaceType SpaceDimension, class T>
class TaylorDerivativesContainer;

/*******************************************************************
 * TaylorDerivativesContainer: Friend Function Declaration
 ******************************************************************/

template<SpaceType SpaceDimension, class T>
bool operator==(const TaylorDerivativesContainer<SpaceDimension,T>& left,
		const TaylorDerivativesContainer<SpaceDimension,T>& right);

template<SpaceType SpaceDimension, class T>
bool operator!= (const TaylorDerivativesContainer<SpaceDimension,T>& left,
		 const TaylorDerivativesContainer<SpaceDimension,T>& right);

template<SpaceType SpaceDimension, class T>
std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<SpaceDimension,T>& Obj);

/*******************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer        *
 ******************************************************/
template<SpaceType SpaceDimension, class T>
class TaylorDerivativesContainer{
 public:
  typedef DerivativeObj<SpaceDimension,T> Derivative;

 private:
  // Variables
  std::vector< Derivative* > DContainer; /* Container of pointers to objects */
 public:
  // Use the provided default Constructor
  TaylorDerivativesContainer(){}
  // Overloaded constructor
  TaylorDerivativesContainer(const unsigned OrderOfReconstruction){ GenerateContainer(OrderOfReconstruction);}
  // Copy constructor
  TaylorDerivativesContainer(const TaylorDerivativesContainer<SpaceDimension,T> & rhs);
  // Destructor
  ~TaylorDerivativesContainer(){ free_memory();}
  // Size()
  unsigned size() const { return DContainer.size();}
  unsigned FirstElem() const { return 0; }
  unsigned LastElem() const { return size()-1;}
  // GenerateContainer()
  void GenerateContainer(const unsigned OrderOfReconstruction); /* Allocate memory for the derivatives and 
								   generate the power combinations */
  // free_memory()
  void free_memory(); 		/* Delete the objects to which the pointers stored in the container point to */
  // IndexOrder(int,int) -> determines the position of the pointer to the element having the powers (p1,p2)
  unsigned IndexOrder(const unsigned p1, const unsigned p2);

  // ComputeSolutionFor( )
  T ComputeSolutionFor(const double DeltaX);
  T ComputeSolutionFor(const double DeltaX, const double DeltaY);
  T ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ);

  /* Overloaded Operators */
  T & operator()(unsigned p1);
  T & operator()(unsigned p1, unsigned p2);
  T & operator()(unsigned p1, unsigned p2, unsigned p3);

  T & operator()(unsigned p1) const;
  T & operator()(unsigned p1, unsigned p2) const;
  T & operator()(unsigned p1, unsigned p2, unsigned p3) const;

  Derivative & operator()(const unsigned position, const bool, const bool, const bool) {return *(DContainer[position]);}
  Derivative & operator()(const unsigned position, const bool, const bool, const bool) const {return *(DContainer[position]);}

  /* Friend functions */
  friend bool operator== <SpaceDimension,T> (const TaylorDerivativesContainer<SpaceDimension,T>& left,
					     const TaylorDerivativesContainer<SpaceDimension,T>& right);

  friend bool operator!= <SpaceDimension,T> (const TaylorDerivativesContainer<SpaceDimension,T>& left,
					     const TaylorDerivativesContainer<SpaceDimension,T>& right);

  friend std::ostream & operator<< <SpaceDimension,T> (std::ostream & os, const TaylorDerivativesContainer<SpaceDimension,T>& Obj);

    // Assignment operator;
  TaylorDerivativesContainer<SpaceDimension,T> & operator=(const TaylorDerivativesContainer<SpaceDimension,T> & rhs);

  // Define the iterator for the DerivativeContainer class
  class iterator;
  friend class iterator;
  class iterator{
    TaylorDerivativesContainer<SpaceDimension,T> & dc;	/* reference to the Container object it iterates thru */
    int index;			/* current location */
  public:
    /* Constructors */
    // set position at the beginning of the container -> index = 0
    iterator(TaylorDerivativesContainer<SpaceDimension,T>& objc) : dc(objc), index(0) { }
    // set position at the end of the container -> index = DContainer.size() - 1;
    iterator(TaylorDerivativesContainer<SpaceDimension,T>& objc, bool) : dc(objc){
      index=dc.DContainer.size(); }

    // Assignment operator
    iterator& operator=(const iterator& rv){
      if( this == &rv) return *this;
      dc = rv.dc;
      index = rv.index;
      return *this;
    }

    // Advance the iterator
    /* Move Forward */
    iterator& operator++( );
    iterator& operator++(int){ //Postfix
      return operator++(); // Use prefix version
    }
    /* Move Backward */
    iterator& operator--( );
    iterator& operator--(int ){ // Postfix
      return operator--(); // Use prefix version
    }
    /* Jump an iterator forward */
    iterator& operator+=(int amount);
    /* Jump an iterator backward */
    iterator& operator-=(int amount);
    /* Create a new iterator that's moved forward */
    iterator operator+(int amount) const;
    /* Create a new iterator that's moved backward */
    iterator operator-(int amount) const;
    /* Comparison */
    bool operator==(const iterator& rv) const;
    bool operator!=(const iterator& rv) const;
    bool operator>=(const iterator& rv) const;
    bool operator<=(const iterator& rv) const;
    bool operator>(const iterator& rv) const;
    bool operator<(const iterator& rv) const;

    /* Access the object from the container */
    DerivativeObj<SpaceDimension,T>* operator->( ) const {
#ifndef __No_Checking__
      require(dc.DContainer[index] != 0, "Zero value "
	      "returned by TaylorDerivativesIterator::operator->()");
#endif
      return dc.DContainer[index];
    }

    /* Friend functions */
    friend std::ostream& operator<<(std::ostream& os, const iterator& it){
      return os << it.index;
    }
  };

  /* Position inside of the container */
  iterator begin() { return iterator(*this); }
  /* Create the "end sentinel" */
  iterator end() { return iterator(*this,true);}
  /* The last element */
  iterator back() { return iterator(*this,true) - 1; }

};

/*******************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer       *
 * Implementation of the Member Functions              *
 ******************************************************/

 /* Copy constructor  */
template<SpaceType SpaceDimension, class T>
inline TaylorDerivativesContainer<SpaceDimension,T>::TaylorDerivativesContainer(const TaylorDerivativesContainer<SpaceDimension,T> & rhs)
{
  // reserve memory for the new container
  DContainer.reserve(rhs.size() );

  for (int i=0; i<rhs.size() ; ++i){
    // allocate memory
    DContainer.push_back(new DerivativeObj<SpaceDimension,T> );
    // copy the value from the LHS
    *(DContainer[i]) = *(rhs.DContainer[i]);
  }
}

/*IndexOrder(unsigned,unsigned) -> determines the position of the pointer to the element having the powers (p1,p2) */
template<SpaceType SpaceDimension, class T>
inline unsigned TaylorDerivativesContainer<SpaceDimension,T>::IndexOrder(const unsigned p1, const unsigned p2)
{

  /* Obs. DContainer.back()->P1() is equal to the order of the reconstruction */
  int sum = 0;

  for (int n=2; n<=p1; ++n){
    sum += (n - 1);
  }

  return (DContainer.back()->P1() + 1)*p1  - sum + p2;
}

/* Allocate memory for the derivatives and power coefficients */
template<SpaceType SpaceDimension, class T>
inline void TaylorDerivativesContainer<SpaceDimension,T>::GenerateContainer(const unsigned OrderOfReconstruction)
{

  // Check if DContainer is empty
  if (DContainer.size() != 0){
    // delete de objects
    free_memory();
    // clear the array
    DContainer.clear();
  }

  int NumberOfDerivatives;	/* The number of derivatives */
  switch(SpaceDimension){
  case OneD:
    /* Determine the number of derivatives */
    NumberOfDerivatives = OrderOfReconstruction+1;
    break;
  case TwoD:
    NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)/2;
    break;
  case ThreeD:
    NumberOfDerivatives = (OrderOfReconstruction+1)*(OrderOfReconstruction+2)*(OrderOfReconstruction+3)/6;
  }

  /* Reserve memory */
  DContainer.reserve(NumberOfDerivatives);
  /* Create the objects */
  for (int i=0; i<NumberOfDerivatives; ++i)
    DContainer.push_back(new DerivativeObj<SpaceDimension,T>(0.0));

  typename std::vector< DerivativeObj<SpaceDimension,T>* >::iterator Iter = DContainer.begin();

  /* First possibility */
#ifdef __Use_Iterator__

  switch(SpaceDimension){
  case OneD:
    for (int p1=0; p1<=OrderOfReconstruction; ++p1){
      (*Iter)->SetPowers(p1);
      Iter++;
    }
    break;
  case TwoD:
    for (int p1=0; p1<=OrderOfReconstruction; ++p1)
      for (int p2=0; p2<=OrderOfReconstruction-p1; ++p2){
	(*Iter)->SetPowers(p1,p2);
	Iter++;
      }
    break;
  case ThreeD:
    for (int p1=0; p1<=OrderOfReconstruction; ++p1)
      for (int p2=0; p2<=OrderOfReconstruction-p1; ++p2)
	for (int p3=0; p3<=OrderOfReconstruction-p1-p2; ++p3){
	  (*Iter)->SetPowers(p1,p2,p3);
	  Iter++;
	}
    break;
  }

#else

  /* Second possibility */
  switch(SpaceDimension){
  case OneD:
    for (int p1=0; p1<=OrderOfReconstruction; ++p1){
      DContainer[p1]->SetPowers(p1);
    }
    break;
  case TwoD:
    for (int p1=0, position=0; p1<=OrderOfReconstruction; ++p1)
      for (int p2=0; p2<=OrderOfReconstruction-p1; ++p2){
	DContainer[position]->SetPowers(p1,p2);
	++position;
      }
    break;
  case ThreeD:
    for (int p1=0, position=0; p1<=OrderOfReconstruction; ++p1)
      for (int p2=0; p2<=OrderOfReconstruction-p1; ++p2)
	for (int p3=0; p3<=OrderOfReconstruction-p1-p2; ++p3){
	  DContainer[position]->SetPowers(p1,p2,p3);
	  ++position;
	}
    break;
  }

#endif
}

/* free_memory() */
template<SpaceType SpaceDimension, class T>
inline void TaylorDerivativesContainer<SpaceDimension,T>::free_memory() {

  /* First possibility */
#ifdef __Use_Iterator__
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::iterator Iter;
  for (Iter=DContainer.begin(); Iter!=DContainer.end(); Iter++){
    delete (*Iter);
    *Iter = NULL;
  }

#else
  /* Second possibility */

  for (int i=0; i<size(); ++i){
    delete DContainer[i];
    DContainer[i] = NULL;
  }

#endif

}

/* Assignment operator = */
template<SpaceType SpaceDimension, class T>
inline TaylorDerivativesContainer<SpaceDimension,T> & 
  TaylorDerivativesContainer<SpaceDimension,T>::operator=(const TaylorDerivativesContainer<SpaceDimension,T> & rhs){

  // !!! If the LHS container already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & rhs) return *this;
  // Check if DContainer is empty
  if (DContainer.size() != 0){
    // delete de objects
    free_memory();
    // clear the array
    DContainer.clear();
  }

  int size = rhs.DContainer.size();
  DContainer.reserve(size);
  for (int i=0; i<size; ++i){
    DContainer.push_back(new DerivativeObj<SpaceDimension,T>);
    *(DContainer[i]) = *(rhs.DContainer[i]);
  }
  return *this;
}

/* Access data */
// operator(unsigned) -> returns the derivative for which p1_ = p1
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()(unsigned p1){
  // This operator is valid only for 1D
#ifndef __No_Checking__
  require(SpaceDimension == OneD, "TaylorDerivativesContainer::operator(unsigned) failed due to difference in the dimensions.\n");
#endif

#ifdef __Use_Iterator__
  /* First possibility */
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;
  
  Iter = DContainer.end();
  Iter--;
  int max_p1 = (*Iter)->P1();
  require(p1 <= max_p1, "TaylorDerivativesContainer::operator(unsigned) failed. Index out of bounds.\n");
  
  Iter = DContainer.begin();
  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1))
      break;
    Iter++;
  }
  return (*Iter)->D();

#else

  /* Second possibility */
  return DContainer[p1]->D();

#endif

}

// operator(unsigned) -> returns the derivative for which p1_ = p1
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()(unsigned p1) const{
  // This operator is valid only for 1D
#ifndef __No_Checking__
  require(SpaceDimension == OneD, "TaylorDerivativesContainer::operator(unsigned) failed due to difference in the dimensions.\n");
#endif

#ifdef __Use_Iterator__ 
  /* First possibility */
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;

  Iter = DContainer.end();
  Iter--;
  int max_p1 = (*Iter)->P1();
  require(p1 <= max_p1, "TaylorDerivativesContainer::operator(unsigned) failed. Index out of bounds.\n");

  Iter = DContainer.begin();
  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1))
      break;
    Iter++;
  }
  return (*Iter)->D();

#else

  /* Second possibility */
  return DContainer[p1]->D();

#endif
}

// operator(unsigned,unsigned) -> returns the derivative for which p1_ = p1 and p2_ = p2
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()(unsigned p1, unsigned p2){
  // This operator is valid only for 2D
#ifndef __No_Checking__
  require(SpaceDimension == TwoD, "TaylorDerivativesContainer::operator(unsigned,unsigned) failed due to difference in the dimensions.\n");
#endif

#ifdef __Use_Iterator__
  /* First possibility */
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;

  Iter = DContainer.end();
  Iter--;
  int max_p1p2 = (*Iter)->P1();
  require( (p1+p2) <= max_p1p2, "TaylorDerivativesContainer::operator(unsigned,unsigned) failed. Index out of bounds.\n");

  Iter = DContainer.begin();
  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1,p2))
      break;
    Iter++;
  }
  return (*Iter)->D();

#else

  /* Second possibility */
  return DContainer[IndexOrder(p1,p2)]->D();

#endif
}

// operator(unsigned,unsigned) -> returns the derivative for which p1_ = p1 and p2_ = p2
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()(unsigned p1, unsigned p2) const{
  // This operator is valid only for 2D
#ifndef __No_Checking__
  require(SpaceDimension == TwoD, "TaylorDerivativesContainer::operator(unsigned,unsigned) failed due to difference in the dimensions.\n");
#endif

#ifdef __Use_Iterator__
  /* First possibility */
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;

  Iter = DContainer.end();
  Iter--;
  int max_p1p2 = (*Iter)->P1();
  require( (p1+p2) <= max_p1p2, "TaylorDerivativesContainer::operator(unsigned,unsigned) failed. Index out of bounds.\n");

  Iter = DContainer.begin();
  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1,p2))
      break;
    Iter++;
  }
  return (*Iter)->D();

#else

  /* Second possibility */
  return DContainer[IndexOrder(p1,p2)]->D();

#endif

}

// operator(unsigned,unsigned,unsigned) -> returns the derivative for which p1_ = p1, p2_ = p2 and p3_ = p3
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()( unsigned p1, unsigned p2, unsigned p3){
  // This operator is valid only for 3D
#ifndef __No_Checking__
  require(SpaceDimension == ThreeD, "TaylorDerivativesContainer::operator(unsigned,unsigned,unsigned) failed due to difference in the dimensions.\n");
#endif
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;

  Iter = DContainer.end();
  Iter--;
  int max_p1p2p3 = (*Iter)->P1();
#ifndef __No_Checking__
  require( (p1+p2+p3) <= max_p1p2p3, "TaylorDerivativesContainer::operator(unsigned,unsigned,unsigned) failed. Index out of bounds.\n");
#endif

  Iter = DContainer.begin();

  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1,p2,p3))
      break;
    Iter++;
  }
  return (*Iter)->D();
}

// operator(unsigned,unsigned,unsigned) -> returns the derivative for which p1_ = p1, p2_ = p2 and p3_ = p3
template<SpaceType SpaceDimension, class T>
inline T & TaylorDerivativesContainer<SpaceDimension,T>::operator()( unsigned p1, unsigned p2, unsigned p3) const{
  // This operator is valid only for 3D
#ifndef __No_Checking__
  require(SpaceDimension == ThreeD, "TaylorDerivativesContainer::operator(unsigned,unsigned,unsigned) failed due to difference in the dimensions.\n");
#endif
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter;

  Iter = DContainer.end();
  Iter--;
  int max_p1p2p3 = (*Iter)->P1();
#ifndef __No_Checking__
  require( (p1+p2+p3) <= max_p1p2p3, "TaylorDerivativesContainer::operator(unsigned,unsigned,unsigned) failed. Index out of bounds.\n");
#endif

  Iter = DContainer.begin();

  while (Iter != DContainer.end()){
    if ((*Iter)->IsPowerEqualTo(p1,p2,p3))
      break;
    Iter++;
  }
  return (*Iter)->D();
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance DeltaX
template<SpaceType SpaceDimension, class T>
inline T TaylorDerivativesContainer<SpaceDimension,T>::ComputeSolutionFor(const double DeltaX){

#ifndef __No_Checking__
  require(SpaceDimension == OneD, "TaylorDerivativesContainer::ComputeSolutionFor( 1D ) failed. Difference in dimensions.\n");
#endif

  // initialize the solution state
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter = DContainer.begin();
  T Solution = (*Iter)->D();
  double value;

  // compute the solution
  for (Iter = DContainer.begin()+1; Iter!=DContainer.end(); Iter++){
    value = std::pow(DeltaX,(double)(*Iter)->P1());
    Solution += value * (*Iter)->D();
  }
  return Solution;
}


// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance DeltaX, DeltaY
template<SpaceType SpaceDimension, class T>
inline T TaylorDerivativesContainer<SpaceDimension,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY){

#ifndef __No_Checking__
  require(SpaceDimension == TwoD, "TaylorDerivativesContainer::ComputeSolutionFor( 2D ) failed. Difference in dimensions.\n");
#endif

  // initialize the solution state
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter = DContainer.begin();
  T Solution = (*Iter)->D();
  double value;

  // compute the solution
  for (Iter = DContainer.begin()+1; Iter!=DContainer.end(); Iter++){
    value = std::pow(DeltaX,(double)(*Iter)->P1()) * std::pow(DeltaY,(double)(*Iter)->P2());
    Solution += value * (*Iter)->D();
  }
  return Solution;
}

// ComputeSolutionFor :-> Compute the solution of the Taylor series expansion for a particular distance DeltaX, DeltaY, DeltaZ
template<SpaceType SpaceDimension, class T>
inline T TaylorDerivativesContainer<SpaceDimension,T>::ComputeSolutionFor(const double DeltaX, const double DeltaY, const double DeltaZ){

#ifndef __No_Checking__
  require(SpaceDimension == ThreeD, "TaylorDerivativesContainer::ComputeSolutionFor( 3D ) failed. Difference in dimensions.\n");
#endif

  // initialize the solution state
  typename std::vector< DerivativeObj<SpaceDimension,T>* >::const_iterator Iter = DContainer.begin();
  T Solution = (*Iter)->D();
  double value;

  // compute the solution
  for (Iter = DContainer.begin()+1; Iter!=DContainer.end(); Iter++){
    value = (std::pow(DeltaX,(double)(*Iter)->P1()) * std::pow(DeltaY,(double)(*Iter)->P2()) *
	     std::pow(DeltaZ,(double)(*Iter)->P3()));
    Solution += value * (*Iter)->D();
  }
  return Solution;
}

// Friend functions
template<SpaceType SpaceDimension, class T>
inline std::ostream & operator<< (std::ostream & os, const TaylorDerivativesContainer<SpaceDimension,T>& Obj){

  int size = Obj.DContainer.size();
  for (int i=0; i<size; ++i)
    os << *(Obj.DContainer[i]);
  return os << std::endl;
}

// operator ==
template<SpaceType SpaceDimension, class T>
inline bool operator==(const TaylorDerivativesContainer<SpaceDimension,T>& left,
		       const TaylorDerivativesContainer<SpaceDimension,T>& right){

  int size = right.DContainer.size();
  bool answer=true;
  if ( left.DContainer.size() != size )
    return !answer;
  for (int i=0; i<size; ++i)
    answer = answer && (*(left.DContainer[i]) == *(right.DContainer[i]));
  return answer;
}

// operator !=
template<SpaceType SpaceDimension, class T>
inline bool operator!=(const TaylorDerivativesContainer<SpaceDimension,T>& left,
		const TaylorDerivativesContainer<SpaceDimension,T>& right){

  return !(left == right);
}

/*****************************************************************
 * TEMPLATIZED CLASS: TaylorDerivativesContainer::iterator       *
 * Implementation of the Member Functions                        *
 ****************************************************************/
// CLASS: DerivativeContainer::iterator
// Iterator Operators


template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator & 
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator++( ) { //Prefix
#ifndef __No_Checking__
  require(++index <= dc.DContainer.size(), "TaylorDerivativesContainer::iterator::operator++ "
	  "moves index out of bounds");
#else
  ++index;
#endif
  return *this;
}

template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator&
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator--( ) { //Prefix
#ifndef __No_Checking__
  require(--index >= 0, "TaylorDerivativesContainer::iterator::operator-- "
	  "moves index out of bounds");
#else
  --index;
#endif
  return *this;
}

template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator&
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator+=(int amount) {
#ifndef __No_Checking__
  require(index + amount < dc.DContainer.size() && index + amount >= 0,
	  "TaylorDerivativesContainer::iterator::operator+= attempt to index out of bounds");
#endif
  index += amount;
  return *this;
}

template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator&
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator-=(int amount) {
#ifndef __No_Checking__
  require(index - amount < dc.DContainer.size() && index - amount >= 0,
	  "TaylorDerivativesContainer::iterator::operator-= attempt to index out of bounds");
#endif
  index -= amount;
  return *this;
}

template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator+(int amount) const{
  iterator ret(*this);
  ret += amount;
  return ret;
}

template<SpaceType SpaceDimension, class T>
inline typename TaylorDerivativesContainer<SpaceDimension,T>::iterator
  TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator-(int amount) const{
  iterator ret(*this);
  ret -= amount;
  return ret;
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator==(const iterator& rv) const{
  return (index == rv.index)&&( dc == rv.dc );
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator!=(const iterator& rv) const{
  return (index != rv.index)&&( dc == rv.dc );
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator>=(const iterator& rv) const{
  return (index >= rv.index)&&( dc == rv.dc );
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator<=(const iterator& rv) const{
  return (index <= rv.index)&&( dc == rv.dc );
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator>(const iterator& rv) const{
  return (index > rv.index)&&( dc == rv.dc );
}

template<SpaceType SpaceDimension, class T>
inline bool TaylorDerivativesContainer<SpaceDimension,T>::iterator::operator<(const iterator& rv) const{
  return (index < rv.index)&&( dc == rv.dc );
}

#endif

#ifndef _SUBGRIDMESH_INCLUDED
#define _SUBGRIDMESH_INCLUDED

#include <vector>
#include "Grid/Grid2D/Cell2D.h"
#include "Common/Containers/PointWiseSolution.h"
#include "include/HeaderData.h"
#include "include/TypeDefinition.h"

template <class NodeType=Node2D, SpaceType SpaceDimension = TwoD, class SolutionType=double >
class SubGridMesh;

/************************************************
*     Friend Functions : ComputationalCell       *
************************************************/
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
bool operator==(const SubGridMesh<NodeType,SpaceDimension,SolutionType>& left,
		const SubGridMesh<NodeType,SpaceDimension,SolutionType>& right);

template<class NodeType, SpaceType SpaceDimension, class SolutionType>
bool operator!=(const SubGridMesh<NodeType,SpaceDimension,SolutionType>& left,
		const SubGridMesh<NodeType,SpaceDimension,SolutionType>& right);

template<class NodeType, SpaceType SpaceDimension, class SolutionType>
std::ostream& operator<< (std::ostream& os,
			  const SubGridMesh<NodeType,SpaceDimension,SolutionType>& SG);

 /****************************************************************************
 * TEMPLATIZED CLASS: SubGridMesh                                           *
 *                                                                          *
 * Container for SubGridPoints.                                             *
 * 
 * VARIABLES:                                                               *
 *                                                                          *
 ***************************************************************************/

template <class NodeType, SpaceType SpaceDimension, class SolutionType>
class SubGridMesh{
 public:
  typedef SubGridMesh<NodeType, SpaceDimension, SolutionType> SubGridMeshType;
  typedef PointWiseSolution<NodeType,SolutionType> SubgridSolutionType;

 private:

  PointWiseSolution<NodeType,SolutionType>*** Ptr;
  vector<int> XYZ;
 public:

  static HeaderData VarNames;

  /* Default constructor and destructor */
  SubGridMesh();
  ~SubGridMesh();
  
  // Copy constructor
  SubGridMesh(const SubGridMeshType & rhs);
  // Assignment operator
  SubGridMeshType & operator=(const SubGridMeshType & rhs);

  /* Memory management functions */
  void allocate();
  void deallocate();
  void SetMesh(const vector<int> &XYZ_);
  void SetMesh(const int & Nx);
  void SetMesh(const int & Nx , const int & Ny);
  void SetMesh(const int & Nx , const int & Ny, const int & Nz);

  /* Field access */
  const int & XPoints(void) const { return XYZ[0]; }
  int & XPoints(void){ return XYZ[0]; }
  const int & YPoints(void) const { return XYZ[1]; }
  int & YPoints(void){ return XYZ[1]; }
  const int & ZPoints(void) const { return XYZ[2]; }
  int & ZPoints(void){ return XYZ[2]; }

  /* Operator overloaded */
  PointWiseSolution<NodeType,SolutionType> & operator( )(const int &i);
  PointWiseSolution<NodeType,SolutionType> & operator( )(const int &i, const int &j);
  PointWiseSolution<NodeType,SolutionType> & operator( )(const int &i, const int &j, const int &k);

  /* Output functions */
  void DefineHeader(const HeaderData & NameOfThePrintedVariables){ VarNames = NameOfThePrintedVariables; }
  void OutputGeometryTecplot(std::ofstream &output_file, const bool Title = false);
  void OutputSolutionTecplot(std::ofstream &output_file, const bool Title = false);
  void OutputSolutionTecplotOneZone(std::ofstream &output_file);
  void OutputSolutionTecplotOneZone(std::ofstream &output_file, const int &SubgridNy);
  void OutputSolutionTecplotOneZone(std::ofstream &output_file, const int &SubgridNy, const int &SubgridNz);
  /* Friend functions */
  friend bool operator== <NodeType,SpaceDimension,SolutionType> (const SubGridMeshType & left,
								 const SubGridMeshType & right);

  friend bool operator!= <NodeType,SpaceDimension,SolutionType> (const SubGridMeshType & left,
								 const SubGridMeshType & right);

  friend std::ostream& operator<< <NodeType,SpaceDimension,SolutionType> (std::ostream& os,
									  const SubGridMeshType& SG);

};


/* Definition of member functions */

// Default constructor
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
SubGridMesh<NodeType,SpaceDimension,SolutionType>::SubGridMesh() :Ptr(NULL)
{
  XYZ.reserve(SpaceDimension);
  for (int i=0; i<SpaceDimension; ++i)
    XYZ.push_back(0);
}

//Copy constructor
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
SubGridMesh<NodeType,SpaceDimension,SolutionType>::SubGridMesh(const SubGridMesh<NodeType, SpaceDimension, SolutionType> & rhs){

  // copy the dimensions
  XYZ = rhs.XYZ;
  // allocate memory
  allocate();
  // copy the elements
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<XYZ[0]; ++i)
      Ptr[0][0][i] = rhs.Ptr[0][0][i];
    break;
  case TwoD:
    for (int i=0; i<XYZ[0]; ++i)
      for (int j=0; j<XYZ[1]; ++j)
	Ptr[0][j][i] = rhs.Ptr[0][j][i];
    break;
  case ThreeD:
    for (int i=0; i<XYZ[0]; ++i)
      for (int j=0; j<XYZ[1]; ++j)
	for (int k=0; k<XYZ[2]; ++k)
	  Ptr[k][j][i] = rhs.Ptr[k][j][i];
    break;
  }
}

// Assignment operator
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  SubGridMesh<NodeType, SpaceDimension, SolutionType> & 
 SubGridMesh<NodeType,SpaceDimension,SolutionType>::operator=(const SubGridMesh<NodeType, SpaceDimension, SolutionType> & rhs) { 
  
  if (this == &rhs) return *this;

  // if the object has memory allocated this is going to be deleted
  if (Ptr!=NULL)
    deallocate();
  // copy the dimensions
  XYZ = rhs.XYZ;
  // allocate memory
  allocate();
  // copy the elements
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<XYZ[0]; ++i)
      Ptr[0][0][i] = rhs.Ptr[0][0][i];
    break;
  case TwoD:
    for (int i=0; i<XYZ[0]; ++i)
      for (int j=0; j<XYZ[1]; ++j)
	Ptr[0][j][i] = rhs.Ptr[0][j][i];
    break;
  case ThreeD:
    for (int i=0; i<XYZ[0]; ++i)
      for (int j=0; j<XYZ[1]; ++j)
	for (int k=0; k<XYZ[2]; ++k)
	  Ptr[k][j][i] = rhs.Ptr[k][j][i];
    break;
  }
  return *this;
}

// Destructor
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
SubGridMesh<NodeType,SpaceDimension,SolutionType>::~SubGridMesh( ) {
  if (Ptr!=NULL)
    deallocate();
}

// Allocate function
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::allocate( ) {
  switch(SpaceDimension){
  case OneD:
    Ptr = new PointWiseSolution<NodeType,SolutionType>** ;
    Ptr[0] = new PointWiseSolution<NodeType,SolutionType>*;
    Ptr[0][0] = new PointWiseSolution<NodeType,SolutionType>[XYZ[0]];
    break;
  case TwoD:
    Ptr = new PointWiseSolution<NodeType,SolutionType>** ;
    Ptr[0] = new PointWiseSolution<NodeType,SolutionType>* [XYZ[1]];
    for (int i=0; i<XYZ[1]; ++i)
      Ptr[0][i] = new PointWiseSolution<NodeType,SolutionType>[XYZ[0]];
    break;
  case ThreeD:
    Ptr = new PointWiseSolution<NodeType,SolutionType>** [XYZ[2]];
    for (int i=0; i<XYZ[2]; ++i)
      Ptr[i] = new PointWiseSolution<NodeType,SolutionType>* [XYZ[1]];
    for (int i=0; i<XYZ[2]; ++i)
      for (int j=0; j<XYZ[1]; ++j)
	Ptr[i][j] = new PointWiseSolution<NodeType,SolutionType>[XYZ[0]];
    break;
  default: 
    std::cout << "In SubGridMesh<>::allocate() the number of dimensions are not properly defined.\n";
  }
}

// Deallocate function
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::deallocate( ) {
  if (Ptr!=NULL){
    switch(SpaceDimension){
    case OneD:
      delete [] Ptr[0][0];
      delete Ptr[0];
      delete Ptr;
      break;
    case TwoD:
      for (int i=0; i<XYZ[1]; ++i)
	delete [] Ptr[0][i];
      delete [] Ptr[0];
      delete Ptr;
      break;
    case ThreeD:
      for (int i=0; i<XYZ[2]; ++i)
	for (int j=0; j<XYZ[1]; ++j)
	  delete [] Ptr[i][j];
      for (int i=0; i<XYZ[2]; ++i)
	delete [] Ptr[i];
      delete [] Ptr;
      break;
    default: 
      std::cout << "In SubGridMesh::deallocate() the number of dimensions are not properly defined.\n";
    }
  }

  Ptr = NULL;
}

// SetMesh(const vector<int> &)
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  void SubGridMesh<NodeType,SpaceDimension,SolutionType>::SetMesh(const vector<int> &XYZ_){
  require(XYZ.size()==XYZ_.size(), "SubGridMesh::SetNumberOfPoints failed\n");
  if(Ptr!=NULL)
    deallocate();
  XYZ = XYZ_;
  allocate();
}

// SetMesh(const int &)
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  void SubGridMesh<NodeType,SpaceDimension,SolutionType>::SetMesh(const int & Nx){
  require(XYZ.size()==OneD, "SubGridMesh::SetNumberOfPoints failed\n");
  if(Ptr!=NULL)
    deallocate();
  XYZ[0] = Nx;
  allocate();
}

// SetMesh(const int &, const int &)
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  void SubGridMesh<NodeType,SpaceDimension,SolutionType>::SetMesh(const int & Nx , const int & Ny){
  require(XYZ.size()==TwoD, "SubGridMesh::SetNumberOfPoints failed\n");
  if(Ptr!=NULL)
    deallocate();
  XYZ[0] = Nx; XYZ[1] = Ny;
  allocate();
}

// SetMesh(const int &, const int &, const int &)
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  void SubGridMesh<NodeType,SpaceDimension,SolutionType>::SetMesh(const int & Nx , const int & Ny, const int & Nz){
  require(XYZ.size()==ThreeD, "SubGridMesh::SetNumberOfPoints failed\n");
  if(Ptr!=NULL)
    deallocate();
  XYZ[0] = Nx; XYZ[1] = Ny; XYZ[2] = Nz;
  allocate();
}

/* Output functions */
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::
  OutputGeometryTecplot(std::ofstream &output_file, const bool Title)
{

  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = SubGridMesh (Node Locations)\n ";
    switch(SpaceDimension){
    case OneD:
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I=" << XYZ[0] << "\n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE) \n";

      break;
    case TwoD:
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \n"
		  << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE DOUBLE) \n";
      break;
    case ThreeD:
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \n"
		  << "\"z\" \n"
		  << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "K = " << XYZ[2] << " \\ \n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE DOUBLE DOUBLE) \n";
      break;
    }
  } else {
    switch(SpaceDimension){
    case OneD:
      output_file << "I=" << XYZ[0] << "\n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE) \n";
      break;
    case TwoD:
      output_file << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE DOUBLE) \n";
      break;
    case ThreeD:
      output_file << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "K = " << XYZ[2] << " \\ \n"
		  << "F = POINT \n"
		  << "DT = (DOUBLE DOUBLE DOUBLE) \n";
      break;
    }
  }

  switch(SpaceDimension){
  case OneD:
    for(int i=0; i<XYZ[0]; ++i)
      output_file << Ptr[0][0][i].GetNode() << "\n";
    break;
  case TwoD:
    for(int j=0; j<XYZ[1]; ++j)
      for(int i=0; i<XYZ[0]; ++i)
	output_file << Ptr[0][j][i].GetNode() << "\n";
    break;
  case ThreeD:
    for(int k=0; k<XYZ[2]; ++k)
      for(int j=0; j<XYZ[1]; ++j)
	for(int i=0; i<XYZ[0]; ++i)
	  output_file << Ptr[k][j][i].GetNode() << "\n";
  }

  output_file << setprecision(6);
}

template<class NodeType, SpaceType SpaceDimension, class SolutionType>
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::
  OutputSolutionTecplot(std::ofstream &output_file, const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = SubGridMesh (Solution At Node Locations)\n ";
    switch(SpaceDimension){
    case OneD:
      output_file << "VARIABLES = \"x\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      output_file << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I=" << XYZ[0] << "\n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
      break;
    case TwoD:
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      output_file << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
      break;
    case ThreeD:
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \n"
		  << "\"z\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      output_file << "ZONE T = \" SubgridMesh \" \\ \n"
		  << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "K = " << XYZ[2] << " \\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
      break;
    }
  } else {
    switch(SpaceDimension){
    case OneD:
      output_file << "I=" << XYZ[0] << "\n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
      break;
    case TwoD:
      output_file << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
      break;
    case ThreeD:
      output_file << "I = " << XYZ[0] << " \\ \n"
		  << "J = " << XYZ[1] << " \\ \n"
		  << "K = " << XYZ[2] << " \\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
      break;
    }
  }

  switch(SpaceDimension){
  case OneD:
    for(int i=0; i<XYZ[0]; ++i)
      output_file << Ptr[0][0][i].GetNode() << "\t" << Ptr[0][0][i].GetValue() << "\n";
    break;
  case TwoD:
    for(int j=0; j<XYZ[1]; ++j)
      for(int i=0; i<XYZ[0]; ++i)
	output_file << Ptr[0][j][i].GetNode() <<  "\t" <<  Ptr[0][j][i].GetValue() << "\n";
    break;
  case ThreeD:
    for(int k=0; k<XYZ[2]; ++k)
      for(int j=0; j<XYZ[1]; ++j)
	for(int i=0; i<XYZ[0]; ++i)
	  output_file << Ptr[k][j][i].GetNode() <<  "\t" << Ptr[k][j][i].GetValue() << "\n";
  }

  output_file << setprecision(6);
}

template<class NodeType, SpaceType SpaceDimension, class SolutionType> inline
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::
  OutputSolutionTecplotOneZone(std::ofstream &output_file)
{
  output_file << setprecision(14);

  for(int i=0; i<XYZ[0]; ++i)
    output_file << Ptr[0][0][i].GetNode() << "\t" << Ptr[0][0][i].GetValue() << "\n";
  output_file << setprecision(6);
}

template<class NodeType, SpaceType SpaceDimension, class SolutionType> inline
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::
  OutputSolutionTecplotOneZone(std::ofstream &output_file, const int &SubgridNy)
{
  output_file << setprecision(14);

  for(int i=0; i<XYZ[0]; ++i)
    output_file << Ptr[0][SubgridNy][i].GetNode() <<  "\t" <<  Ptr[0][SubgridNy][i].GetValue() << "\n";

  output_file << setprecision(6);
}

template<class NodeType, SpaceType SpaceDimension, class SolutionType> inline
void SubGridMesh<NodeType,SpaceDimension,SolutionType>::
  OutputSolutionTecplotOneZone(std::ofstream &output_file, const int &SubgridNy, const int &SubgridNz)
{
  output_file << setprecision(14);

  for(int i=0; i<XYZ[0]; ++i)
    output_file << Ptr[SubgridNz][SubgridNy][i].GetNode() <<  "\t" <<  Ptr[SubgridNz][SubgridNy][i].GetValue() << "\n";

  output_file << setprecision(6);
}

/* Overloaded operators */
// bracket ( int )
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  PointWiseSolution<NodeType,SolutionType> & 
  SubGridMesh<NodeType,SpaceDimension,SolutionType>::operator( )(const int &i){
  require(SpaceDimension == OneD, "SubGridMesh::operator( ) failed." 
	  " You are not working in 1D. Wrong number of parameters passed!.\n");
  require((i<XYZ[0])&&(i>=0),"SubGridMesh::operator( ) failed." 
	  " The indices are out of bounds! \n");
  return Ptr[0][0][i];
}

// bracket ( int, int )
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  PointWiseSolution<NodeType,SolutionType> & 
  SubGridMesh<NodeType,SpaceDimension,SolutionType>::operator( )(const int &i,const int &j){
  require(SpaceDimension == TwoD, "SubGridMesh::operator( ) failed." 
	  " You are not working in 2D. Wrong number of parameters passed!.\n");
  require((i<XYZ[0])&&(i>=0)&&(j<XYZ[1])&&(j>=0), "SubGridMesh::operator( ) failed."
	  " The indices are out of bounds! \n"); 

  return Ptr[0][j][i];
}

// bracket ( int, int , int)
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
  PointWiseSolution<NodeType,SolutionType> & 
  SubGridMesh<NodeType,SpaceDimension,SolutionType>::operator( )(const int &i,const int &j,const int &k){
  require(SpaceDimension == ThreeD, "SubGridMesh::operator( ) failed." 
	  " You are not working in 3D. Wrong number of parameters passed!.\n");
  return Ptr[k][j][i];
}

/* Friend functions */
// operator ==
template<class NodeType, SpaceType SpaceDimension, class SolutionType> inline
bool operator==(const SubGridMesh<NodeType,SpaceDimension,SolutionType>& left,
		const SubGridMesh<NodeType,SpaceDimension,SolutionType>& right)
{
  bool result = true;
  result = result && (left.XYZ==right.XYZ);
  if (result){
    switch(SpaceDimension){
    case OneD:
      for (int i=0; i<left.XYZ[0]; ++i)
	result = result && (left.Ptr[0][0][i] == right.Ptr[0][0][i]);
      break;
    case TwoD:
      for (int i=0; i<left.XYZ[0]; ++i)
	for (int j=0; j<left.XYZ[1]; ++j)
	   result = result && (left.Ptr[0][j][i] == right.Ptr[0][j][i]);
      break;
     case ThreeD:
       for (int i=0; i<left.XYZ[0]; ++i)
	 for (int j=0; j<left.XYZ[1]; ++j)
	   for (int k=0; k<left.XYZ[2]; ++k)
	     result = result && (left.Ptr[k][j][i] == right.Ptr[k][j][i]);
       
       break;
    }
  }
  return result;
}

// operator !=
template<class NodeType, SpaceType SpaceDimension, class SolutionType> inline
bool operator!=(const SubGridMesh<NodeType,SpaceDimension,SolutionType>& left,
		const SubGridMesh<NodeType,SpaceDimension,SolutionType>& right)
{
  return !(left == right);
}

// operator <<
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
std::ostream& operator<< (std::ostream& os,
			  const SubGridMesh<NodeType,SpaceDimension,SolutionType>& SG)
{
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<SG.XYZ[0]; ++i)
      os << SG.Ptr[0][0][i];
    break;
  case TwoD:
    for (int i=0; i<SG.XYZ[0]; ++i)
      for (int j=0; j<SG.XYZ[1]; ++j)
	os << SG.Ptr[0][j][i];
    break;
  case ThreeD:
    for (int i=0; i<SG.XYZ[0]; ++i)
      for (int j=0; j<SG.XYZ[1]; ++j)
	for (int k=0; k<SG.XYZ[2]; ++k)
	  os << SG.Ptr[k][j][i];
    break;
  }
  return os;
}

// Static variables
template<class NodeType, SpaceType SpaceDimension, class SolutionType>
HeaderData SubGridMesh<NodeType,SpaceDimension,SolutionType>::VarNames("Solution");

#endif

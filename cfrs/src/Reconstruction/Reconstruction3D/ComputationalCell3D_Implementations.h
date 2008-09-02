// ComputationalCell3D_Implementations.h defines the implementation of member functions of template ComputationalCell3D

/*******************************************************
 * TEMPLATIZED CLASS: ComputationalCell                *
 * Implementation of member functions                  *
 ******************************************************/

// Constructor
template< class GeometryType,class SolutionType> inline
ComputationalCell<ThreeD,GeometryType,SolutionType>::ComputationalCell():
  geom(),U_cell(0.0),TD(),SubgridSolution(),rings(0),RO(0),
     FinalOrder(0),ErrorRecL1Norm(0.0),ErrorRecL2Norm(0.0),ErrorRecMaxNorm(0.0){
  NbSubgridPoints.reserve(ThreeD);
}

// Copy Constructor
template< class GeometryType,class SolutionType> inline
ComputationalCell<ThreeD,GeometryType,SolutionType>::ComputationalCell( const ComputationalCell & rhs):
  geom(rhs.geom), U_cell(rhs.U_cell),rings(rhs.rings),RO(rhs.RO),FinalOrder(rhs.FinalOrder),
     ErrorRecL1Norm(rhs.ErrorRecL1Norm), ErrorRecL2Norm(rhs.ErrorRecL2Norm), ErrorRecMaxNorm(rhs.ErrorRecMaxNorm),
     NbSubgridPoints(rhs.NbSubgridPoints),SubgridSolution(rhs.SubgridSolution){ }

// Generate cell in 3D;
template< class GeometryType,class SolutionType>
void ComputationalCell<ThreeD,GeometryType,SolutionType>::GenerateCell(const int NbX, const int NbY, const int NbZ, const int OrderOfRec){

  RO = OrderOfRec;
  NbSubgridPoints.reserve(ThreeD);
  NbSubgridPoints.push_back(NbX);
  NbSubgridPoints.push_back(NbY);
  NbSubgridPoints.push_back(NbZ);
  // generate the container of Taylor derivatives
  TD.GenerateContainer(OrderOfRec);
  // generate SubgridMesh
  SubgridSolution.SetMesh(NbSubgridPoints[0],NbSubgridPoints[1],NbSubgridPoints[2]);
}

// Set Rings
template< class GeometryType, class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::SetRings(void){
  int ND;
  ND = TD.size();

  // ThreeD
  switch(ND){
  case 1:  rings = 0; break;
  case 4:  rings = 1; break;
  case 10: rings = 1; break;
  case 20: rings = 2; break;
  case 35: rings = 2; break;
  default:
    rings = -1;
  }
}

// SetCell Parameters
template< class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int & RO_){
  require(ThreeD <= SubGridPoints.size(), "ComputationalCell::SetCell() failed. The size of the provided vector is too small for the required dimension\n");
  // copy geometry
  SetCell(geom_);

  // generate SubgridMesh and TaylorDerivatives container
  GenerateCell(SubGridPoints[0], SubGridPoints[1], SubGridPoints[2], RO_);

  // compute the geometric coefficients
  ComputeGeometricCoefficients();

  // Set number of rings
  SetRings();

  // Set geometry of the subgrid
  SetSubGridGeometry();
}

template< class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  SetCell (const std::vector<int> &SubGridPoints, const int & RO_){
  require(ThreeD <= SubGridPoints.size(), "ComputationalCell::SetCell() failed. The size of the provided vector is too small for the required dimension\n");

  // generate SubgridMesh and TaylorDerivatives container
  GenerateCell(SubGridPoints[0], SubGridPoints[1], SubGridPoints[2], RO_);

  // Set number of rings
  SetRings();
}

// SetSubGridGeometry()
template< class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::SetSubGridGeometry(){
  Vector XStart;
  
  std::cout << "Not implemented yet\n";

}


// ComputeGeometricCoefficients()
/******************************************************************
 * Computes the geometric coefficients by computing integral      *
 * over the volume of the cell of the polynomial function having  *
 * the form                                                       *
 *       ((x-xCC)^n)*((y-yCC)^m) and divide by aria A             * 
******************************************************************/
template<class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>
::ComputeGeometricCoefficients( ){

  int l,m,n;
  double DummyParam;

  /* Create the polynomial function for the cell */
  GeneralizedPolynomialFunctionOfThreeVariables Polynom(0,0,0,geom.xc.x,geom.xc.y,geom.xc.z);

//#ifdef __Use_Iterator__
//  /* Create iterator */
//  GeomCoeffIterator Iter = GeomCoeff;
//
//  for (Iter = GeomCoeff.begin(); Iter!=GeomCoeff.end(); Iter++){
//    /* Compute the geometric coefficient */
//    l = Iter->P1();
//    m  = Iter->P2();
//    n = Iter->P3();
//    Polynom.ChangePowersTo(l,m,n);
//    Iter->D() = IntegrateOverTheCell(Polynom,14,DummyParam)/geom.V();
//  }
//#else


  /* When the method of integration is used the moments associated with the first order are not quite accurate */
  for (int p1 = GeomCoeff.FirstElem(); p1<=GeomCoeff.LastElem(); ++p1){
    /* Compute the geometric coefficient */
    Polynom.ChangePowersTo(GeomCoeff(p1,true,true,true).P1(),GeomCoeff(p1,true,true,true).P2(),GeomCoeff(p1,true,true,true).P3());
    GeomCoeff(p1,true,true,true).D() = IntegrateOverTheCell(Polynom,14,DummyParam)/geom.V();
  }
}



// UpdateSubgridSolution()
template< class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::UpdateSubgridSolution(){

  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j)
      for (int k=0; k<NbSubgridPoints[2]; ++k){
	SubgridSolution(i,j,k).SetValue( SolutionAt(SubgridSolution(i,j,k).GetNode()) );
      }
}

// template<typename FO>
// UpdateExactSubgridSolution(const FO)
template<class GeometryType,class SolutionType>
  template<typename FO> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  UpdateExactSubgridSolution(const FO FuncObj)
{
  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j)
      for (int k=0; k<NbSubgridPoints[2]; ++k){
	SubgridSolution(i,j,k).SetValue( ExactSolutionAt(SubgridSolution(i,j,k).GetNode(), FuncObj) );
      }
}


//SolutionAtCoordinates(double & , double & , double & )
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<ThreeD,GeometryType,SolutionType>::
SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord){

  double DifferenceX, DifferenceY, DifferenceZ;
  DifferenceX = X_Coord - geom.xc.x;
  DifferenceY = Y_Coord - geom.xc.y;
  DifferenceZ = Z_Coord - geom.xc.z;

  return TD.ComputeSolutionFor(DifferenceX,DifferenceY,DifferenceZ);
}

template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<ThreeD,GeometryType,SolutionType>::
SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord, int VarPosition){
  return SolutionAtCoordinates(X_Coord,Y_Coord,Z_Coord)[VarPosition];
}


// Compute solution : SolutionAt()
template< class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<ThreeD,GeometryType,SolutionType>::
  SolutionAt(double & deltaX, double & deltaY, double & deltaZ){
  return TD.ComputeSolutionFor(deltaX,deltaY,deltaZ);
}

template< class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<ThreeD,GeometryType,SolutionType>::
  SolutionAt(Node & node){

  Vector difference;
  
  difference = node.X - geom.xc;
  return  TD.ComputeSolutionFor(difference.x,difference.y,difference.z);
}

// template<typename FunctionObject>
// void ComputeReconstructionError(const FO FuncObj)
/*******************************************************************
Computes the error of the reconstructed function in all three norms
********************************************************************/
template<class GeometryType,class SolutionType>
  template<typename FO>
void ComputationalCell<ThreeD,GeometryType,SolutionType>::ComputeReconstructionError(const FO FuncObj){

  /* integrate the error for the L1 norm */
  ErrorRecL1Norm = IntegrateOverTheCell(error_function(FuncObj,
						       wrapped_member_function(this,
									       &SolutionAtCoordinates,
									       ErrorRecL1Norm)),
					10, ErrorRecL1Norm);

  /* integrate the error for the L2 norm*/
  ErrorRecL2Norm = IntegrateOverTheCell(square_error_function(FuncObj,
							      wrapped_member_function(this,
										      &SolutionAtCoordinates,
										      ErrorRecL2Norm)),
					10, ErrorRecL2Norm);

  /* generate 101 points over the cell domain and chose the maximum error */
  ErrorRecMaxNorm = SolutionType(0.0);
  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j)
      for (int k=0; k<NbSubgridPoints[2]; ++k){
	ErrorRecMaxNorm = max(ErrorRecMaxNorm,
			      error_function(FuncObj,
					     wrapped_member_function(this,
								     &SolutionAtCoordinates,
								     ErrorRecL1Norm))(SubgridSolution(i,j,k).GetNode().x(),
										      SubgridSolution(i,j,k).GetNode().y(),
										      SubgridSolution(i,j,k).GetNode().z()));
      }
}

// OutputSubgridTecplot()
/* Does the output of the coordinates of the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
		       const int & kCell, const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell3D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \n"
		<< "\"z\" \n"
		<< "ZONE T = \" Cell [" << iCell <<"," << jCell <<"," << kCell << "] \" \\ \n";
  } else {
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell <<"," << kCell << "] \" \\ \n";
  }

  SubgridSolution.OutputGeometryTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionTecplot()
/* Does the output of the solution at the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			const int & kCell, const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell3D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \n"
		<< "\"z\" \n";
    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell <<"," << kCell << "] \" \\ \n";

  } else {
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell <<"," << kCell << "] \" \\ \n";
  }

  SubgridSolution.OutputSolutionTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionCellTecplot()
/* Does the output of the first derivative solution at the cell center */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			    const int & kCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell3D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n"
		<< "\"z\" \n";
    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell <<"," << kCell  << "] \" \\ \n";
  }

  output_file << geom.xc << "\t" << TD(0,0,0) << "\n";

  output_file << setprecision(6);
}

// OutputSolutionCellTecplotOneZone()
/* Does the output of the first derivative solution at the cell center */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<ThreeD,GeometryType,SolutionType>::
  OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int &SubgridNy, const int &SubgridNz)
{
  SubgridSolution.OutputSolutionTecplotOneZone(output_file,SubgridNy,SubgridNz);
}

// Friend functions
// Operator ==
template< class GeometryType,class SolutionType> inline
bool operator==(const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
		const ComputationalCell<ThreeD,GeometryType,SolutionType>& right){

  return (left.geom==right.geom)&&(left.U_cell==right.U_cell)&&(left.TD==right.TD)&&(left.NbSubgridPoints==right.NbSubgridPoints);
}

// Operator !=
template< class GeometryType,class SolutionType> inline
bool operator!=(const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
		const ComputationalCell<ThreeD,GeometryType,SolutionType>& right){
  return !(left==right);
}

// Operator <<
template< class GeometryType,class SolutionType> inline
std::ostream& operator<< (std::ostream& os, const ComputationalCell<ThreeD,GeometryType,SolutionType>& Obj){
  os << std::endl << "geom=" << Obj.geom << std::endl ;
  os << "Average Solution=" << Obj.U_cell << std::endl;
  os << " Derivatives:" << std::endl << Obj.TD;
  os << "rings=" << Obj.rings << std::endl;
  os << "ReconstructionOrder=" << Obj.RO << std::endl;
  os << "FinalOrder=" << Obj.FinalOrder << std::endl;
  os << "ReconstructionError(L1, L2, LMax)=" << Obj.ErrorRecL1Norm 
     << "\t" << Obj.ErrorRecL2Norm << "\t" << Obj.ErrorRecMaxNorm << std::endl;
  os << "#SubgridPoints:" ;
  for (int i=0; i<Obj.NbSubgridPoints.size(); ++i)
    os << "\tdir " << i+1 << "-> " << Obj.NbSubgridPoints[i] << " points,";
  os << std::endl;
  os << "SubgridSolution:" << std::endl << Obj.SubgridSolution << std::endl;
}

// Static variables
template<class GeometryType,class SolutionType>
HeaderData ComputationalCell<ThreeD,GeometryType,SolutionType>::VarNames("Solution");

/* Include full specialization of template member functions */
#include "ComputationalCell3D_Specializations.h"

// ComputationalCell3D.h defines the partial specialization of ComputationalCell template for the
// SpaceDimension = ThreeD
// Parameters' definition is given in ComputationalCellTemplate.h

#ifndef _COMPUTATIONALCELL3D_INCLUDED
#define _COMPUTATIONALCELL3D_INCLUDED

// Define the 3D ComputationalCell class : Partial specialization of ComputationalCell
template< class GeometryType, class SolutionType>
class ComputationalCell<ThreeD,GeometryType,SolutionType>;

/************************************************
*     Friend Functions : ComputationalCell      *
************************************************/
template<class GeometryType, class SolutionType>
bool operator==(const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
		const ComputationalCell<ThreeD,GeometryType,SolutionType>& right);

template<class GeometryType, class SolutionType>
bool operator!=(const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
		const ComputationalCell<ThreeD,GeometryType,SolutionType>& right);

template<class GeometryType, class SolutionType>
std::ostream& operator<< (std::ostream& os, const ComputationalCell<ThreeD,GeometryType,SolutionType>& Obj);

 /****************************************************************************
 * TEMPLATIZED CLASS: ComputationalCell 3D                                   *
 * Partial Specialization for ThreeD                                         *
 ****************************************************************************/
template < class GeometryType, class SolutionType>
class ComputationalCell<ThreeD, GeometryType, SolutionType>{
 public:

  typedef typename GeometryTraits<ThreeD>::NodeType Node;
  typedef typename GeometryTraits<ThreeD>::VectorType Vector;
  typedef TaylorDerivativesContainer<ThreeD,SolutionType> DerivativesContainer;
  typedef typename DerivativesContainer::iterator DerivIterator;
  typedef TaylorDerivativesContainer<ThreeD,double> GeometricIntegrals;
  typedef typename GeometricIntegrals::iterator GeomCoeffIterator;
  typedef SubGridMesh<Node,ThreeD,SolutionType> SubGridType;
  typedef ComputationalCell<ThreeD, GeometryType, SolutionType> CompCellType;

  static const int NumberOfVariables = SolutionParameters<SolutionType>::NUM_OF_VARIABLES;

 private:

  /* Cell Variables */
  GeometryType geom;
  SolutionType U_cell;
  SolutionType MCC;		/* Multiple-Correlation Coefficient */
  DerivativesContainer TD;
  GeometricIntegrals GeomCoeff;
  int rings;
  int RO;
  int FinalOrder;
  SolutionType ErrorRecL1Norm;
  SolutionType ErrorRecL2Norm;
  SolutionType ErrorRecMaxNorm;  

  /* SubDomain Variables */
  vector<int> NbSubgridPoints;
  SubGridType SubgridSolution;

  void AdjustNumberOfSubgridPoints(void){};

 public:

  /* Static Variable */
  static HeaderData VarNames;

  // Constructor
  ComputationalCell();
  // Copy Constructor
  ComputationalCell( const ComputationalCell & rhs);
  // Generate cell
  void GenerateCell (const int NbX, const int NbY, const int OrderOfRec);
  void GenerateCell (const int NbX, const int NbY, const int NbZ, const int OrderOfRec);
  void SetRings (void);
  void SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int &RO_);
  void SetCell (const std::vector<int> &SubGridPoints, const int & RO_);
  void SetCell (const GeometryType & geom_){ geom = geom_;}

  // Field access
  const Vector & CellCenter( ) {return geom.xc;}
  GeometryType& CellGeometry( ) { return geom;}
  SolutionType & CellSolution( ) {return U_cell;}
  DerivativesContainer & CellDeriv() {return TD;}
  DerivIterator BeginD() {return TD.begin();}
  DerivIterator EndD() {return TD.end();}
  DerivIterator BackD() {return TD.back();}
  SolutionType & CellDeriv(int i, int j, int k) {return TD(i,j,k);}
  const int & CellRings() const {return rings;}
  const int & CellRecOrder() const {return RO;}
  int & CellFOrder() {return FinalOrder;}

  const SolutionType & CellErrorL1() const { return ErrorRecL1Norm;}
  SolutionType & CellErrorL1() { return ErrorRecL1Norm;}
  const SolutionType & CellErrorL2() const { return ErrorRecL2Norm;}
  SolutionType & CellErrorL2() { return ErrorRecL2Norm;}
  const SolutionType & CellErrorMax() const { return ErrorRecMaxNorm;}
  SolutionType & CellErrorMax() { return ErrorRecMaxNorm;}

  // Added by RR:
  const GeometricIntegrals & CellGeomCoeff() const { return GeomCoeff; }
  GeometricIntegrals & CellGeomCoeff() {return GeomCoeff;}
  const double & CellGeomCoeff(const int & p1, const int & p2, cont int & p3) {return GeomCoeff(p1,p2,p3);}
#ifndef __Use_Iterator__
  const double & CellGeomCoeff(const int position, const bool, const bool, const bool) {return GeomCoeff(position,true,true,true).D();}
#endif
  //


  const int & iSubgridPoints() const {return NbSubgridPoints[0];}
  const int & jSubgridPoints() const {return NbSubgridPoints[1];}
  const int & kSubgridPoints() const {return NbSubgridPoints[2];}

  void SetSubGridGeometry(void);
  void UpdateSubgridSolution(void);
  SolutionType SolutionAt(double & deltaX, double & deltaY, double & deltaZ);
  SolutionType SolutionAt(Node & node);
  SolutionType SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord){
    std::cout<< "CompCell3D::SolutionAtCoordinates not defined!";
  }
  // Compute exact solution
  template<typename FO>
    SolutionType ExactSolutionAt(Node & node, const FO FuncObj){
    return FuncObj(node.x(),node.y(),node.z());
  }
  template<typename FO>
    void UpdateExactSubgridSolution(const FO FuncObj);

  /* Operating functions */
  template<typename PointerToTestFunction>
    void SetInitialSolution(PointerToTestFunction PtrFunc){}

  // ComputeReconstructionError()
  template<typename FO>
    void ComputeReconstructionError(const FO FuncObj);

  // IntegrateOverTheCell()
  template<typename FO, class ReturnType>
    ReturnType IntegrateOverTheCell(const FO FuncObj, const int & digits, ReturnType _dummy_param){
    return _dummy_param ;
  }

  void OutputMeshNodesTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			      const int & kCell,const bool Title=false){ };
  void OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			    const int & kCell,const bool Title=false);
  void OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			     const int & kCell,const bool Title=false);
  void OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
				 const int & kCell,const bool Title=false);
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int &SubgridNy, 
					const int &SubgridNz);
  void OutputPWC(std::ofstream &output_file, const int &SubgridNy, const int &SubgridNz){ };
  void OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const int &SubgridNy, 
					       const int &SubgridNz);

  void OutputReconstructedSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
					  const int & kCell,const bool Title=false){};

  void OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file, int _dummy1_, int _dummy2_){ };

  // Friend functions
  // == operator
  friend bool operator== <GeometryType,SolutionType>
    (const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
     const ComputationalCell<ThreeD,GeometryType,SolutionType>& right);

  // != operator
  friend bool operator!= <GeometryType,SolutionType>
    (const ComputationalCell<ThreeD,GeometryType,SolutionType>& left,
     const ComputationalCell<ThreeD,GeometryType,SolutionType>& right);

  // << operator
  friend std::ostream& operator<< <GeometryType,SolutionType>
    (std::ostream& os, const ComputationalCell<ThreeD,GeometryType,SolutionType>& Obj);
};

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

#endif

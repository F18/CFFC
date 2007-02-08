// ComputationalCell1D.h defines the partial specialization of ComputationalCell template for the
// SpaceDimension = OneD
// Parameters' definition is given in ComputationalCellTemplate.h

#ifndef _COMPUTATIONALCELL1D_INCLUDED
#define _COMPUTATIONALCELL1D_INCLUDED

// Define the 1D ComputationalCell class : Partial specialization of ComputationalCell
template< class GeometryType, class SolutionType>
class ComputationalCell<OneD,GeometryType,SolutionType>;

/************************************************
*     Friend Functions : ComputationalCell1D     *
************************************************/
template<class GeometryType, class SolutionType>
bool operator==(const ComputationalCell<OneD,GeometryType,SolutionType>& left,
		const ComputationalCell<OneD,GeometryType,SolutionType>& right);

template<class GeometryType, class SolutionType>
bool operator!=(const ComputationalCell<OneD,GeometryType,SolutionType>& left,
		const ComputationalCell<OneD,GeometryType,SolutionType>& right);

template<class GeometryType, class SolutionType>
std::ostream& operator<< (std::ostream& os, const ComputationalCell<OneD,GeometryType,SolutionType>& Obj);

/*******************************************************
 * CLASS Template: ComputationalCell                   *
 * Partial Specialization for OneD                     *
 ******************************************************/
template <class GeometryType, class SolutionType>
class ComputationalCell<OneD, GeometryType, SolutionType>{
 public:
  
  typedef typename GeometryTraits<OneD>::NodeType Node;
  typedef typename GeometryTraits<OneD>::VectorType Vector;
  typedef TaylorDerivativesContainer<OneD,SolutionType> DerivativesContainer;
  typedef TaylorDerivativesContainer<OneD,double> GeometricIntegrals;
#ifdef __Use_Iterator__
  typedef typename DerivativesContainer::iterator DerivIterator;
  typedef typename GeometricIntegrals::iterator GeomCoeffIterator;
#else
  typedef typename DerivativesContainer::Derivative  Derivative;
#endif
  typedef SubGridMesh<Node,OneD,SolutionType> SubGridType;
  typedef typename SubGridType::SubgridSolutionType SubGridSolutionType;
  typedef ComputationalCell<OneD, GeometryType, SolutionType> CompCellType;
  
  static const int NumberOfVariables = SolutionParameters<SolutionType>::NUM_OF_VARIABLES;

 private:

  /* Cell Variables */
  GeometryType geom;
  SolutionType U_cell;
  SolutionType MCC;		/* Multiple-Correlation Coefficient */
  DerivativesContainer TD;
  DerivativesContainer TD_WENO;
  DerivativesContainer TD_FirstOrder; /* container of derivatives for first order reconstruction */
  GeometricIntegrals GeomCoeff;
  int rings;
  int RO;
  int FinalOrder;
  SolutionType ErrorRecL1Norm;
  SolutionType ErrorRecL2Norm;
  SolutionType ErrorRecMaxNorm;

  /* Monotonicity Information */
  int LimitedCell;

  /* SubDomain Variables */
  vector<int> NbSubgridPoints;
  SubGridType SubgridSolution;

  void AdjustNumberOfSubgridPoints(void);

 public:

  /* Static Variable */
  static HeaderData VarNames;

  // Constructor
  ComputationalCell();
  // Copy Constructor
  ComputationalCell( const ComputationalCell & rhs);
  // Generate cell
  void GenerateCell (const int NbX, const int OrderOfRec);
  void SetRings (void);
  void SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int & RO_);
  void SetCell (const std::vector<int> & SubGridPoints, const int & RO_){};
  void SetCell (const GeometryType & geom_) { geom = geom_;}

  // Field access

  const double & CellCenter(void) const {return geom.x;}
  const double & CellDelta (void) {return geom.dx;}
  const double & CellDomain(void) {return geom.dx;}
  const GeometryType& CellGeometry(void) { return geom;}

  const SolutionType& CellSolution(void) const {return U_cell;}
  SolutionType & CellSolution(void) {return U_cell;}
  const double & CellSolution(int & VarPosition) const;
  double & CellSolution(int & VarPosition);

  const SolutionType& CellMCC(void) const {return MCC;}
  SolutionType & CellMCC(void) {return MCC;}
  const double & CellMCC(int & VarPosition) const;
  double & CellMCC(int & VarPosition);

  const DerivativesContainer & CellDeriv(void) const {return TD;}
  DerivativesContainer & CellDeriv(void) {return TD;}
  const SolutionType & CellDeriv(int i) const {return TD(i);}
  SolutionType & CellDeriv(int i) {return TD(i);}

  const DerivativesContainer & CellDerivWENO(void) const {return TD_WENO;}
  DerivativesContainer & CellDerivWENO(void) {return TD_WENO;}
  const SolutionType & CellDerivWENO(int i) const {return TD_WENO(i);}
  SolutionType & CellDerivWENO(int i) {return TD_WENO(i);}

  const DerivativesContainer & CellDerivFirstOrder(void) const {return TD_FirstOrder;}
  DerivativesContainer & CellDerivFirstOrder(void) {return TD_FirstOrder;}
  const SolutionType & CellDerivFirstOrder(int i) const {return TD_FirstOrder(i);}
  SolutionType & CellDerivFirstOrder(int i) {return TD_FirstOrder(i);}

  const SolutionType & CellLimiter(void) const {return TD.Limiter();}
  SolutionType & CellLimiter(void) {return TD.Limiter();}

  const int & UnfitReconstructionFlag(void) const {return LimitedCell; }
  int & UnfitReconstructionFlag(void){return LimitedCell; }

#ifndef __Use_Iterator__
  Derivative & CellDeriv(const unsigned position, const bool, const bool, const bool) {return TD(position,true,true,true);}
  Derivative & CellDerivFirstOrder(const unsigned position,
				   const bool, const bool, const bool) {return TD_FirstOrder(position,true,true,true);}
#endif
  int NumberOfTaylorDerivatives() const {return TD.size();}

#ifdef __Use_Iterator__
  DerivIterator BeginD() {return TD.begin();}
  DerivIterator EndD() {return TD.end();}
  DerivIterator BackD() {return TD.back();}
#endif

  const int & CellRings() {return rings;}
  const int & CellRecOrder() {return RO;}

  const int & CellFOrder() const {return FinalOrder;}
  int & CellFOrder() {return FinalOrder;}

  const SolutionType & CellErrorL1() const { return ErrorRecL1Norm;}
  SolutionType & CellErrorL1() { return ErrorRecL1Norm;}
  const SolutionType & CellErrorL2() const { return ErrorRecL2Norm;}
  SolutionType & CellErrorL2() { return ErrorRecL2Norm;}
  const SolutionType & CellErrorMax() const { return ErrorRecMaxNorm;}
  SolutionType & CellErrorMax() { return ErrorRecMaxNorm;}

  const GeometricIntegrals & CellGeomCoeff() const { return GeomCoeff; }
  GeometricIntegrals & CellGeomCoeff() {return GeomCoeff;}
  const double & CellGeomCoeff(const unsigned & p1) {return GeomCoeff(p1);}

  const SubGridSolutionType & CellSubgridSolution(const int &i) const {return SubgridSolution(i);}
  SubGridSolutionType & CellSubgridSolution(const int &i) {return SubgridSolution(i);}

#ifndef __Use_Iterator__
  const double & CellGeomCoeff(const int position, const bool, const bool, const bool) {
    return GeomCoeff(position,true,true,true).D();
  }
#endif

  unsigned iSubgridPoints() const {return NbSubgridPoints[0];}
  unsigned jSubgridPoints() const {return 0;}
  unsigned kSubgridPoints() const {return 0;}

  // Operating functions
  void SetSubGridGeometry(void);
  void ComputeGeometricCoefficients(void);
  void UpdateSubgridSolution(void);
  CompCellType Reflexion(double & ReflexionCentre){ return ComputationalCell(this);};
  SolutionType SolutionAt(const double & deltaX);
  SolutionType SolutionAtCoordinates(double & X_Coord);
  // Compute exact solution
  template<typename FO>
    SolutionType ExactSolutionAt(Node & node, const FO FuncObj){
    return FuncObj(node);
  }
  template<typename FO>
    void UpdateExactSubgridSolution(const FO FuncObj);

  void OutputMeshNodesTecplot(std::ofstream &output_file,const int & iCell,const bool Title=false){ };

  void OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const bool Title=false);
  void OutputSubgridTecplot(std::ofstream &output_file,const int & iCell, const int & _dummy1_,
			    const int & _dummy2_, const bool Title=false){
    OutputSubgridTecplot(output_file,iCell,Title);
  };
  void OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const bool Title=false);
  void OutputSolutionTecplot(std::ofstream &output_file,const int & iCell, const int & _dummy1_,
			     const int & _dummy2_, const bool Title=false){
    OutputSolutionTecplot(output_file,iCell,Title);
  };
  void OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const bool Title=false);
  void OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell, const int & _dummy1_,
				 const int & _dummy2_, const bool Title=false){
    OutputSolutionCellTecplot(output_file,iCell,Title);
  };
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file);
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int & _dummy1_, const int & _dummy2_){
    OutputSolutionCellTecplotOneZone(output_file);
  };
  void OutputPWC(std::ofstream &output_file);
  void OutputPWC(std::ofstream &output_file, const int & _dummy1_, const int & _dummy2_){
    OutputPWC(output_file);
  };
  void OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const int & _dummy1_, const int & _dummy2_){  };
  void OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file, int _dummy1_, int _dummy2_){ };


  void OutputReconstructedSolutionTecplot(std::ofstream &output_file,const int & iCell,const bool Title=false);
  void OutputReconstructedSolutionTecplot(std::ofstream &output_file, const int & iCell, const int & _dummy1_,
					  const int & _dummy2_, const bool Title=false){
    OutputReconstructedSolutionTecplot(output_file,iCell,Title);
  };

  template<typename PointerToTestFunction>
    void SetInitialSolution(PointerToTestFunction PtrFunc);

  // IntegrateOverTheCell()
  template<typename FO, class ReturnType>
    ReturnType IntegrateOverTheCell(const FO FuncObj, const int & digits, ReturnType _dummy_param){
    return AdaptiveGaussianQuadrature(FuncObj,geom.x-0.5*geom.dx,geom.x+0.5*geom.dx,digits,_dummy_param);
  }

  // ComputeReconstructionError()
  template<typename FO>
    void ComputeReconstructionError(const FO FuncObj);

  // Friend functions
  // == operator
  friend bool operator== <GeometryType,SolutionType> (const ComputationalCell<OneD,GeometryType,SolutionType>& left,
						      const ComputationalCell<OneD,GeometryType,SolutionType>& right);

  // != operator
  friend bool operator!= <GeometryType,SolutionType> (const ComputationalCell<OneD,GeometryType,SolutionType>& left,
						      const ComputationalCell<OneD,GeometryType,SolutionType>& right);

  // << operator
  friend std::ostream& operator<< <GeometryType,SolutionType>
    (std::ostream& os, const ComputationalCell<OneD,GeometryType,SolutionType>& Obj);
};

/*******************************************************
 * CLASS Template: ComputationalCell                   *
 * Partial Specialization for OneD                     *
 * Implementation of member functions                  *
 ******************************************************/

// Constructor
template<class GeometryType,class SolutionType> inline
ComputationalCell<OneD, GeometryType,SolutionType>::ComputationalCell():
  geom(),U_cell(0.0),TD(),TD_WENO(),TD_FirstOrder(),SubgridSolution(),rings(0),RO(0),
     FinalOrder(0),ErrorRecL1Norm(0.0),ErrorRecL2Norm(0.0),ErrorRecMaxNorm(0.0), LimitedCell(OFF){
  NbSubgridPoints.reserve(OneD);
}

// Copy Constructor
template<class GeometryType,class SolutionType> inline
ComputationalCell<OneD,GeometryType,SolutionType>::ComputationalCell( const ComputationalCell & rhs):
  geom(rhs.geom), U_cell(rhs.U_cell),rings(rhs.rings),RO(rhs.RO),FinalOrder(rhs.FinalOrder),
     ErrorRecL1Norm(rhs.ErrorRecL1Norm), ErrorRecL2Norm(rhs.ErrorRecL2Norm), ErrorRecMaxNorm(rhs.ErrorRecMaxNorm),
     NbSubgridPoints(rhs.NbSubgridPoints),SubgridSolution(rhs.SubgridSolution){ }

// AdjustNumberOfSubgridPoints -> this function is called before the subgrid is generated
//                                but after the number of subgrid points is initialized
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::AdjustNumberOfSubgridPoints(void){

  /* Set the proper number of points in the x direction */
  /* X dir */
  if (NbSubgridPoints[0] < 3){
    NbSubgridPoints[0] = 3;
  } else if (NbSubgridPoints[0]%2 != 1){
    NbSubgridPoints[0] ++; //only odd numbers
  }
}

// Generate cell in 1D;
template<class GeometryType,class SolutionType>
void ComputationalCell<OneD,GeometryType,SolutionType>::GenerateCell(const int NbX, const int OrderOfRec){
  RO = OrderOfRec;
  FinalOrder = RO;
  NbSubgridPoints.push_back(NbX);

  // adjust the number of subgrid points
  AdjustNumberOfSubgridPoints();

  // generate SubgridMesh
  SubgridSolution.SetMesh(NbSubgridPoints[0]);

  // generate the container of Taylor derivatives
  TD.GenerateContainer(OrderOfRec);
  TD_WENO.GenerateContainer(OrderOfRec);
  TD_FirstOrder.GenerateContainer(1); /* OrderOfRec = 1 */
  
  // generate the container of the geometric coefficients
  GeomCoeff.GenerateContainer(OrderOfRec);
}

// Set Rings
template<class GeometryType, class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::SetRings(void){
  int ND;
  ND = TD.size();

  switch(ND){
  case 1:
    rings = 0;
    break;
  case 2:
    rings = 1;
    break;
  case 3:
    rings = 2;
    break;
  case 4:
    rings = 2;
    break;
  case 5:
    rings = 3;
    break;
  case 6:
    rings = 3;
    break;
  case 7:
    rings = 4;
    break;
  default:
    rings = -1;
  }
}

// SetCell Parameters
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int & RO_){
  require(OneD <= SubGridPoints.size(), "ComputationalCell::SetCell() failed. The size of the provided vector is too small for the required dimension\n");
  // copy geometry
  SetCell(geom_);

  // generate SubgridMesh and TaylorDerivatives container
  GenerateCell(SubGridPoints[0], RO_);

  // compute the geometric coefficients
  ComputeGeometricCoefficients();

  // Set number of rings
  SetRings();

  // Set geometry of the subgrid
  SetSubGridGeometry();
}

// SetSubGridGeometry()
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::SetSubGridGeometry(){
  double XStart, deltaSG;

  deltaSG = geom.dx/(NbSubgridPoints[0]-1);
  XStart = geom.x - geom.dx/2;
  for (int i=0; i<=NbSubgridPoints[0]-1; ++i){
    SubgridSolution(i).SetNode(XStart + i*deltaSG);
  }
}

// ComputeGeometricCoefficients()
/******************************************************************
 * Computes the geometric coefficients by computing integral      *
 * over the volume of the cell of the polynomial function having  *
 * the form                                                       *
 *       ((x-xCC)^n) and divide by length dx                      * 
******************************************************************/
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>
::ComputeGeometricCoefficients(void){

  int n;

  /* Create the polynomial function for the cell */
  GeneralizedPolynomialFunctionOfOneVariable Polynom(0,geom.x);

#ifdef __Use_Iterator__
  double DummyParam = 0.0;
  /* Create iterator */
  GeomCoeffIterator Iter = GeomCoeff;

  for (Iter = GeomCoeff.begin(); Iter!=GeomCoeff.end(); Iter++){
    /* Compute the geometric coefficient */
    n = Iter->P1();
    Polynom.ChangePowerTo(n);
    Iter->D() = IntegrateOverTheCell(Polynom,14,DummyParam)/geom.dx;
  }
#else
  for (int p1=GeomCoeff.FirstElem(); p1<=GeomCoeff.LastElem(); ++p1){
    Polynom.ChangePowerTo(p1);
    GeomCoeff(p1) = IntegrateOverTheCell(Polynom,14,GeomCoeff(p1))/geom.dx;
  }
#endif
}

// UpdateSubgridSolution()
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::UpdateSubgridSolution(){

  double DistanceToPoint;
  for (int i=0; i<NbSubgridPoints[0]; ++i){
    DistanceToPoint = SubgridSolution(i).GetNode() - geom.x;
    SubgridSolution(i).SetValue(TD.ComputeSolutionFor(DistanceToPoint));
  }
}

// template<typename FO>
// UpdateExactSubgridSolution(const FO)
template<class GeometryType,class SolutionType>
  template<typename FO> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  UpdateExactSubgridSolution(const FO FuncObj)
{
  for (int i=0; i<NbSubgridPoints[0]; ++i){
    SubgridSolution(i).SetValue( ExactSolutionAt(SubgridSolution(i).GetNode(), FuncObj) );
  }
}

// Compute solution : SolutionAt()
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<OneD,GeometryType,SolutionType>::
  SolutionAt(const double & deltaX){
  return TD.ComputeSolutionFor(deltaX);
}

//SolutionAtCoordinates(double &)
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<OneD,GeometryType,SolutionType>::
  SolutionAtCoordinates(double & X_Coord){

  double DifferenceX;
  DifferenceX = X_Coord - geom.x;

  return TD.ComputeSolutionFor(DifferenceX);
}

// CellSolution(int & VarPosition) const;
template<class GeometryType,class SolutionType> inline
const double & ComputationalCell<OneD,GeometryType,SolutionType>::
  CellSolution(int & VarPosition) const {

  return U_cell[VarPosition];
}

// CellSolution(int & VarPosition);
template<class GeometryType,class SolutionType> inline
double & ComputationalCell<OneD,GeometryType,SolutionType>::
  CellSolution(int & VarPosition) {
  
  return U_cell[VarPosition];
}

// CellMCC(int & VarPosition) const;
template<class GeometryType,class SolutionType> inline
const double & ComputationalCell<OneD,GeometryType,SolutionType>::
  CellMCC(int & VarPosition) const {

  return MCC[VarPosition];
}

// CellMCC(int & VarPosition);
template<class GeometryType,class SolutionType> inline
double & ComputationalCell<OneD,GeometryType,SolutionType>::
  CellMCC(int & VarPosition) {
  
  return MCC[VarPosition];
}

/* // Reflexion(double & ReflexionCentre); */
/* template<class GeometryType,class SolutionType> inline */
/* CompCellType ComputationalCell<OneD,GeometryType,SolutionType>:: */
/*   CompCellType Reflexion(double & ReflexionCentre){   */
  
/*   return ; */
/* } */

// template<typename PointerToTestFunction>
// SetInitialSolution(PointerToTestFunction)
/***************************************************************
Sets the value of the U_cell based on the integration of PtrFunc
over the domain of the cell
*****************************************************************/
template<class GeometryType,class SolutionType>
  template<typename PointerToTestFunction> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::SetInitialSolution(PointerToTestFunction PtrFunc)
{
  /* compute the average solution */
  U_cell = IntegrateOverTheCell(PtrFunc,6,U_cell)/geom.dx;
  /* update the first coefficient of the Taylor expansion */
  TD(0) = U_cell;
  TD_WENO(0) = U_cell;
  TD_FirstOrder(0) = U_cell;
}

// template<typename FunctionObject>
// void ComputeReconstructionError(const FO FuncObj)
/*******************************************************************
Computes the error of the reconstructed function in all three norms
********************************************************************/
template<class GeometryType,class SolutionType>
  template<typename FO>
void ComputationalCell<OneD,GeometryType,SolutionType>::ComputeReconstructionError(const FO FuncObj){

  /* construct the error function for the L1 norm*/
  _Error_<FO,CompCellType, SolutionType> ErrorFunction(FuncObj,this);
  /* integrate the error */
  ErrorRecL1Norm = IntegrateOverTheCell(ErrorFunction,10,ErrorRecL1Norm);

  /* construct the error function for the L2 norm*/
  _Error_Square_<FO,CompCellType, SolutionType> ErrorFunctionSquare(FuncObj,this);
  /* integrate the error */
  ErrorRecL2Norm = IntegrateOverTheCell(ErrorFunctionSquare,10,ErrorRecL2Norm);

  /* generate 101 points over the cell domain and chose the maximum error */
  ErrorRecMaxNorm = 0.0;
  double DeltaMaxNorm = CellDelta()/100;
  double StartPoint = CellCenter() - 0.5*CellDelta();
  for (int i=0; i<=100; ++i){
    ErrorRecMaxNorm = max(ErrorRecMaxNorm,ErrorFunction(StartPoint+i*DeltaMaxNorm));
  }

}

// OutputSubgridTecplot()
/* Does the output of the coordinates of the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell1D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "ZONE T = \" Cell [" << iCell << "] \" \\ \n";
  } else {
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n";
  }

  SubgridSolution.OutputGeometryTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionTecplot()
/* Does the output of the solution at the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell1D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n";

    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n";

  } else {
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n";
  }

  SubgridSolution.OutputSolutionTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionCellTecplot()
/* Does the output of the first derivative solution at the cell center */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell1D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n";

    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n";
    VarNames.PrintDataTypeTecplot(output_file,OneD);
  }

  output_file << geom.x << "\t" << TD(0) << "\n";

  output_file << setprecision(6);
}

// OutputSolutionCellTecplotOneZone()
/* Does the output of the whole domain solution at the nodes in only one zone */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputSolutionCellTecplotOneZone(std::ofstream &output_file)
{
  SubgridSolution.OutputSolutionTecplotOneZone(output_file);
}

// OutputPWC()
/* Outputs the piecewise constant solution at the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputPWC(std::ofstream &output_file)
{

  output_file << setprecision(14);

  for(int i=0; i<iSubgridPoints(); ++i)
    output_file << SubgridSolution(i).GetNode() <<  "\t" <<  CellSolution() << "\n";

  output_file << setprecision(6);
}

// OutputReconstructedSolutionTecplot()
/* Outputs the reconstructed functions */
template<class GeometryType,class SolutionType>
void ComputationalCell<OneD,GeometryType,SolutionType>::
  OutputReconstructedSolutionTecplot(std::ofstream &output_file,const int & iCell,const bool Title)
{

  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell1D (Reconstructed Solution)\n ";
    output_file << "VARIABLES = \"x\" \\ \n";
    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n"
		<< "I=" << (1 + 2*CellRings())*iSubgridPoints() << "\n"	/* the total number of subgrid points
									   in the stencil used for reconstruction */
		<< "F = POINT \n";
    VarNames.PrintDataTypeTecplot(output_file,OneD);

  } else {
    output_file << "ZONE T = \" Cell [" << iCell << "] \" \\ \n"
      		<< "I=" << (1 + 2*CellRings())*iSubgridPoints() << "\n" 
		<< "F = POINT \n";
    VarNames.PrintDataTypeTecplot(output_file,OneD);
  }

  int cell, SubgridPoint;

  for(cell=iCell-CellRings(); cell<=iCell+CellRings(); ++cell){
    for(SubgridPoint=0; SubgridPoint<iSubgridPoints(); ++SubgridPoint){
      output_file << SubgridSolution(SubgridPoint).GetNode() << "\n";
    }
  }

  output_file << setprecision(6);
}

// Friend functions
// Operator ==
template<class GeometryType,class SolutionType> inline
bool operator==(const ComputationalCell<OneD,GeometryType,SolutionType>& left,
		const ComputationalCell<OneD,GeometryType,SolutionType>& right){

  return (left.geom==right.geom)&&(left.U_cell==right.U_cell)&&(left.TD==right.TD)&&(left.NbSubgridPoints==right.NbSubgridPoints);
}

// Operator !=
template<class GeometryType,class SolutionType> inline
bool operator!=(const ComputationalCell<OneD,GeometryType,SolutionType>& left,
		const ComputationalCell<OneD,GeometryType,SolutionType>& right){
  return !(left==right);
}

// Operator <<
template<class GeometryType,class SolutionType> inline
std::ostream& operator<< (std::ostream& os, const ComputationalCell<OneD,GeometryType,SolutionType>& Obj){
  os << std::endl << "geom=" << Obj.geom << std::endl ;
  os << "Average Solution=" << Obj.U_cell << std::endl;
  os << "Derivatives:" << std::endl << Obj.TD;
  os << "Geometric Coefficients:" << std::endl << Obj.GeomCoeff;
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

  return os;
}

// Static variables
template<class GeometryType,class SolutionType>
HeaderData ComputationalCell<OneD,GeometryType,SolutionType>::VarNames("Solution");

// Specialization for SolutionType == double
template<> inline
const double & ComputationalCell<OneD,Cell1D_NonUniform,double>::
 CellSolution(int &) const {
  return U_cell;
}

template<> inline
double & ComputationalCell<OneD,Cell1D_NonUniform,double>::
 CellSolution(int &){
  return U_cell;
}


template<> inline
const double & ComputationalCell<OneD,Cell1D_Uniform,double>::
 CellSolution(int &) const {
  return U_cell;
}

template<> inline
double & ComputationalCell<OneD,Cell1D_Uniform,double>::
 CellSolution(int &){
  return U_cell;
}

#endif

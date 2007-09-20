// ComputationalCell2D.h defines the partial specialization of ComputationalCell template for the
// SpaceDimension = TwoD
// Parameters' definition is given in ComputationalCellTemplate.h

#ifndef _COMPUTATIONALCELL2D_INCLUDED
#define _COMPUTATIONALCELL2D_INCLUDED

// Define the 2D ComputationalCell class : Partial specialization of ComputationalCell
template<class GeometryType, class SolutionType>
class ComputationalCell<TwoD,GeometryType,SolutionType>;

/************************************************
*     Friend Functions : ComputationalCell       *
************************************************/
template< class GeometryType, class SolutionType>
bool operator==(const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
		const ComputationalCell<TwoD,GeometryType,SolutionType>& right);

template< class GeometryType, class SolutionType>
bool operator!=(const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
		const ComputationalCell<TwoD,GeometryType,SolutionType>& right);

template< class GeometryType, class SolutionType>
std::ostream& operator<< (std::ostream& os, const ComputationalCell<TwoD,GeometryType,SolutionType>& Obj);

 /****************************************************************************
 * TEMPLATIZED CLASS: ComputationalCell2D                                   *
 * Partial Specialization for TwoD                                          *
 ***************************************************************************/
template < class GeometryType, class SolutionType>
class ComputationalCell<TwoD, GeometryType, SolutionType>{
 public:

  typedef typename GeometryTraits<TwoD>::NodeType Node;
  typedef typename GeometryTraits<TwoD>::VectorType Vector;
  typedef TaylorDerivativesContainer<TwoD,SolutionType> DerivativesContainer;
  typedef TaylorDerivativesContainer<TwoD,double> GeometricIntegrals;
#ifdef __Use_Iterator__
  typedef typename DerivativesContainer::iterator DerivIterator;
  typedef typename GeometricIntegrals::iterator GeomCoeffIterator;
#else
  typedef typename DerivativesContainer::Derivative  Derivative;
#endif
  typedef SubGridMesh<Node,TwoD,SolutionType> SubGridType;

  typedef ComputationalCell<TwoD, GeometryType, SolutionType> CompCellType;
  typedef SolutionType (CompCellType::*MemberFunctionType2D)(double &, double &);

  static const int NumberOfVariables = SolutionParameters<SolutionType>::NUM_OF_VARIABLES;

 private:

  /* Cell Variables */
  GeometryType geom;
  SolutionType U_cell;
  SolutionType MCC;		/* Multiple-Correlation Coefficient */
  DerivativesContainer TD;
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

  // Construct cell functions
  void GenerateCell (const int NbX, const int NbY, const int OrderOfRec);
  void SetRings (void);
  void SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int &RO_);
  void SetCell (const std::vector<int> &SubGridPoints, const int & RO_);
  void SetCell (const GeometryType & geom_){ geom = geom_;}
 
  // Field access
  //  Vector & CellCenter() {return geom.xc;}
  const Vector & CellCenter( ) const {return geom.xc;}

  const GeometryType& CellGeometry() const {return geom;}
  GeometryType& CellGeometry( ) { return geom;}

  const double CellDomain() const {return geom.A();}

  const SolutionType& CellSolution() const {return U_cell;}
  SolutionType & CellSolution( ) {return U_cell;}
  const double & CellSolution(int & VarPosition) const;
  double & CellSolution(int & VarPosition);

  const SolutionType& CellMCC(void) const {return MCC;}
  SolutionType & CellMCC(void) {return MCC;}
  const double & CellMCC(int & VarPosition) const;
  double & CellMCC(int & VarPosition);

  int NumberOfTaylorDerivatives() const {return TD.size();}

  /* CellDeriv() */
  const DerivativesContainer & CellDeriv() const {return TD;}
  DerivativesContainer & CellDeriv() {return TD;}
  const SolutionType & CellDeriv(int i, int j) const {return TD(i,j);}
  SolutionType & CellDeriv(int i, int j) {return TD(i,j);}

#ifndef __Use_Iterator__
  Derivative & CellDeriv(const unsigned position, const bool, const bool, const bool) { return TD(position,true,true,true);}
  const double & CellDeriv(const int p1, const int p2, const int Variable) const { return TD(p1,p2,Variable);}
  double & CellDeriv(const int p1, const int p2, const int Variable) { return TD(p1,p2,Variable);}

  /* CellDerivFirstOrder() */
  const DerivativesContainer & CellDerivFirstOrder(void) const {return TD_FirstOrder;}
  DerivativesContainer & CellDerivFirstOrder(void) {return TD_FirstOrder;}
  const SolutionType & CellDerivFirstOrder(int i, int j) const {return TD_FirstOrder(i,j);}
  SolutionType & CellDerivFirstOrder(int i, int j) {return TD_FirstOrder(i,j);}

  Derivative & CellDerivFirstOrder(const unsigned position, const bool, const bool, const bool) {
    return TD_FirstOrder(position,true,true,true);
  }
  const double & CellDerivFirstOrder(const int p1, const int p2, const int Variable) const {
    return TD_FirstOrder(p1,p2,Variable);
  }
  double & CellDerivFirstOrder(const int p1, const int p2, const int Variable) {
    return TD_FirstOrder(p1,p2,Variable);
  }
#endif

#ifdef __Use_Iterator__
  DerivIterator BeginD() {return TD.begin();}
  DerivIterator EndD() {return TD.end();}
  DerivIterator BackD() {return TD.back();}
#endif

  /* Monotonicity Information */
  const SolutionType & CellLimiter(void) const {return TD.Limiter();}
  SolutionType & CellLimiter(void) {return TD.Limiter();}

  const int & UnfitReconstructionFlag(void) const {return LimitedCell; }
  int & UnfitReconstructionFlag(void){return LimitedCell; }

  const int & CellRings() const {return rings;}
  const int & CellRecOrder() const {return RO;}

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
  const double & CellGeomCoeff(const int & p1, const int & p2) {return GeomCoeff(p1,p2);}

#ifndef __Use_Iterator__
  const double & CellGeomCoeff(const int position, const bool, const bool, const bool) {return GeomCoeff(position,true,true,true).D();}
#endif

  const int & iSubgridPoints() const {return NbSubgridPoints[0];}
  const int & jSubgridPoints() const {return NbSubgridPoints[1];}
  int kSubgridPoints() const {return 0;}

  // Operating functions
  void SetSubGridGeometry(void);
  void ComputeGeometricCoefficients(void);
  void UpdateSubgridSolution(void);
  SolutionType SolutionAt( double & deltaX, double & deltaY);
  SolutionType SolutionAtCoordinates( double & X_Coord, double & Y_Coord);
  double SolutionAtCoordinates( double & X_Coord, double & Y_Coord, int VarPosition);
  SolutionType SolutionAt(Node & node);
  // Compute exact solution
  template<typename FO>
    SolutionType ExactSolutionAt(Node & node, const FO FuncObj){
    return FuncObj(node.x(),node.y());
  }
  template<typename FO>
    void UpdateExactSubgridSolution(const FO FuncObj);

  void OutputMeshNodesTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title=false);
  void OutputMeshNodesTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			      const int & _dummy_,const bool Title=false){
    OutputMeshNodesTecplot(output_file,iCell,jCell,Title);
  };
  void OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title=false);
  void OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			    const int & _dummy_,const bool Title=false){
    OutputSubgridTecplot(output_file,iCell,jCell,Title);
  };
  void OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title=false);
  void OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
			     const int & _dummy_,const bool Title=false){
    OutputSolutionTecplot(output_file,iCell,jCell,Title);
  };
  void OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,
				 const int & jCell,const bool Title=false);
  void OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
				 const int & _dummy_,const bool Title=false){
    OutputSolutionCellTecplot(output_file,iCell,jCell,Title);
  };
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int &SubgridNy);
  void OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int &SubgridNy, const int & _dummy_){
    OutputSolutionCellTecplotOneZone(output_file,SubgridNy);
  };

  void OutputPWC(std::ofstream &output_file, const int &SubgridNy);
  void OutputPWC(std::ofstream &output_file, const int &SubgridNy, const int & _dummy_){
    OutputPWC(output_file,SubgridNy);
  };

  void OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const int &SubgridNy);
  void OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const int &SubgridNy, const int & _dummy_){
    OutputSmoothnessIndicatorAtNodesTecplot(output_file,SubgridNy);
  };
  void OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file, int Face, int GaussPoint);


  void OutputReconstructedSolutionTecplot(std::ofstream &output_file,
					  const int & iCell,const int & jCell,const bool Title=false){};
  void OutputReconstructedSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,
					  const int & _dummy_,const bool Title=false){
    OutputReconstructedSolutionTecplot(output_file,iCell,jCell,Title);
  };

  template<typename PointerToTestFunction>
    void SetInitialSolution(PointerToTestFunction PtrFunc);

  // ComputeReconstructionError()
  template<typename FO>
    void ComputeReconstructionError(const FO FuncObj);

  // IntegrateOverTheCell()
  template<typename FO, class ReturnType>
    ReturnType IntegrateOverTheCell(FO FuncObj, const int & digits, const ReturnType & _dummy_param){

    Node2D SW(geom.NodeSW()), NW(geom.NodeNW());
    Node2D SE(geom.NodeSE()), NE(geom.NodeNE());

    return QuadrilateralQuadrature(FuncObj, 
				   SW, NW, NE, SE,
   				   digits,_dummy_param);
  }

  // Friend functions
  // == operator
  friend bool operator== <GeometryType,SolutionType>
    (const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
     const ComputationalCell<TwoD,GeometryType,SolutionType>& right);

  // != operator
  friend bool operator!= <GeometryType,SolutionType>
    (const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
     const ComputationalCell<TwoD,GeometryType,SolutionType>& right);

  // << operator
  friend std::ostream& operator<< <GeometryType,SolutionType>
    (std::ostream& os, const ComputationalCell<TwoD,GeometryType,SolutionType>& Obj);
};

/* Include the implementations of template member functions */
#include "ComputationalCell2D_Implementations.h"

#endif

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
  typedef SolutionType (CompCellType::*MemberFunctionType1D)(const double &);
  
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
  const double & CellDeriv(const int p1, const int Variable) const { return TD(p1,Variable);}
  double & CellDeriv(const int p1, const int Variable) { return TD(p1,Variable);}
  Derivative & CellDerivFirstOrder(const unsigned position,
				   const bool, const bool, const bool) {return TD_FirstOrder(position,true,true,true);}
  const double & CellDerivFirstOrder(const int p1, const int Variable) const {
    return TD_FirstOrder(p1,Variable);
  }
  double & CellDerivFirstOrder(const int p1, const int Variable) {
    return TD_FirstOrder(p1,Variable);
  }

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
  SolutionType SolutionAtCoordinates(const double & X_Coord);
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
    return AdaptiveGaussianQuadrature(FuncObj,geom.x-0.5*geom.dx,geom.x+0.5*geom.dx,_dummy_param,digits);
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


/* Include the implementations of template member functions */
#include "ComputationalCell1D_Implementations.h"

#endif

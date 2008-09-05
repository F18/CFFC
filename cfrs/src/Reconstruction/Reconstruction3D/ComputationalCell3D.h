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
  typedef TaylorDerivativesContainer<ThreeD,double> GeometricIntegrals;
  typedef typename DerivativesContainer::Derivative  Derivative;
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
  void ComputeGeometricCoefficients (void);

  // Field access
  const Vector & CellCenter( ) {return geom.xc;}
  GeometryType& CellGeometry( ) { return geom;}
  const double CellDomain() {return geom.V();}
  //  const double CellDomain() const {return geom.V();}

//  DerivIterator BeginD() {return TD.begin();}
//  DerivIterator EndD() {return TD.end();}
//  Deriviterator BackD() {return TD.back();}
  
  const int & CellRings() const {return rings;}
  const int & CellRecOrder() const {return RO;}
  int & CellFOrder() {return FinalOrder;}

  const SolutionType & CellErrorL1() const { return ErrorRecL1Norm;}
  SolutionType & CellErrorL1() { return ErrorRecL1Norm;}
  const SolutionType & CellErrorL2() const { return ErrorRecL2Norm;}
  SolutionType & CellErrorL2() { return ErrorRecL2Norm;}
  const SolutionType & CellErrorMax() const { return ErrorRecMaxNorm;}
  SolutionType & CellErrorMax() { return ErrorRecMaxNorm;}

  const GeometricIntegrals & CellGeomCoeff() const { return GeomCoeff; }
  GeometricIntegrals & CellGeomCoeff() {return GeomCoeff;}
  const double & CellGeomCoeff(const int & p1, const int & p2, const int & p3) {return GeomCoeff(p1,p2,p3);}
#ifndef __Use_Iterator__
  const double & CellGeomCoeff(const int position, const bool, const bool, const bool) {return GeomCoeff(position,true,true,true).D();}
#endif

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
  const SolutionType & CellDeriv(int i, int j, int k) const {return TD(i,j,k);}
  SolutionType & CellDeriv(int i, int j, int k) {return TD(i,j,k);}
  const DerivativesContainer & CellDeriv() const {return TD;}
  DerivativesContainer & CellDeriv() {return TD;}
#ifndef __Use_Iterator__
  Derivative & CellDeriv(const unsigned position) { return TD(position,true,true,true);}
  const double & CellDeriv(const int p1, const int p2, const int p3, const int Variable) const { return TD(p1,p2,p3,Variable);}
  double & CellDeriv(const int p1, const int p2,  const int p3, const int Variable) { return TD(p1,p2,p3,Variable);}

//  /* CellDerivFirstOrder() */
//  const DerivativesContainer & CellDerivFirstOrder(void) const {return TD_FirstOrder;}
//  DerivativesContainer & CellDerivFirstOrder(void) {return TD_FirstOrder;}
//  const SolutionType & CellDerivFirstOrder(int i, int j, int k) const {return TD_FirstOrder(i,j,k);}
//  SolutionType & CellDerivFirstOrder(int i, int j, int k) {return TD_FirstOrder(i,j,k);}
//
//  Derivative & CellDerivFirstOrder(const unsigned position, const bool, const bool, const bool) {
//    return TD_FirstOrder(position,true,true,true);
//  }
//  const double & CellDerivFirstOrder(const int p1, const int p2, const int p3, const int Variable) const {
//    return TD_FirstOrder(p1,p2,p3,Variable);
//  }
//  double & CellDerivFirstOrder(const int p1, const int p2, const int p3, const int Variable) {
//    return TD_FirstOrder(p1,p2,p3,Variable);
//  }
#endif
  //


  const int & iSubgridPoints() const {return NbSubgridPoints[0];}
  const int & jSubgridPoints() const {return NbSubgridPoints[1];}
  const int & kSubgridPoints() const {return NbSubgridPoints[2];}

  void SetSubGridGeometry(void);
  void UpdateSubgridSolution(void);
  SolutionType SolutionAt(double & deltaX, double & deltaY, double & deltaZ);
  SolutionType SolutionAt(Node & node);
  SolutionType SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord);
  SolutionType SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord, int VarPosition);

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
    ReturnType IntegrateOverTheCell(FO FuncObj, const int & digits, const ReturnType & _dummy_param){

    Node3D SWTop(geom.NodeSWTop()), NWTop(geom.NodeNWTop());
    Node3D SETop(geom.NodeSETop()), NETop(geom.NodeNETop());
    Node3D SWBot(geom.NodeSWBot()), NWBot(geom.NodeNWBot());
    Node3D SEBot(geom.NodeSEBot()), NEBot(geom.NodeNEBot());

    return AdaptiveGaussianQuadrature(FuncObj, NWBot.X.x, NEBot.X.x, NWBot.X.y, SWBot.X.y,
                                      NWBot.X.z, NWTop.X.z, digits,_dummy_param);
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

/* Include the implementations of template member functions */
#include "ComputationalCell3D_Implementations.h"

#endif

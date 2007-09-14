// ComputationalCell2D_Specializations.h defines specializations of member functions of template ComputationalCell2D

/*******************************************************
 * CLASS Template: ComputationalCell2D                 *
 * Specialization of member functions                  *
 ******************************************************/

// Specialization for SolutionType == double & GeometryType == Cell2D_Quad/Cell2D_Cartesian

/* CellSolution() for Cell2D_Quad */
template<> inline
const double & ComputationalCell<TwoD,Cell2D_Quad,double>::
 CellSolution(int &) const {
  return U_cell;
}

template<> inline
double & ComputationalCell<TwoD,Cell2D_Quad,double>::
 CellSolution(int &){
  return U_cell;
}

/* CellSolution() for Cell2D_Cartesian */
template<> inline
const double & ComputationalCell<TwoD,Cell2D_Cartesian,double>::
 CellSolution(int &) const {
  return U_cell;
}

template<> inline
double & ComputationalCell<TwoD,Cell2D_Cartesian,double>::
 CellSolution(int &){
  return U_cell;
}

/* SolutionAtCoordinates() for Cell2D_Quad */
template<> inline
double ComputationalCell<TwoD,Cell2D_Quad,double>::
  SolutionAtCoordinates(double & X_Coord, double & Y_Coord, int VarPosition){
  return SolutionAtCoordinates(X_Coord, Y_Coord);
}

template<> inline
double ComputationalCell<TwoD,Cell2D_Cartesian,double>::
  SolutionAtCoordinates(double & X_Coord, double & Y_Coord, int VarPosition){
  return SolutionAtCoordinates(X_Coord, Y_Coord);
}

/* CellMCC() for Cell2D_Quad */
template<> inline
double & ComputationalCell<TwoD,Cell2D_Quad,double>::
 CellMCC(int &){
  return MCC;
}

template<> inline
const double & ComputationalCell<TwoD,Cell2D_Quad,double>::
 CellMCC(int &) const{
  return MCC;
}

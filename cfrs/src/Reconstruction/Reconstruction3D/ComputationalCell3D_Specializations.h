// ComputationalCell2D_Specializations.h defines specializations of member functions of template ComputationalCell2D

/*******************************************************
 * CLASS Template: ComputationalCell2D                 *
 * Specialization of member functions                  *
 ******************************************************/

// Specialization for SolutionType == double & GeometryType == Cell3D_Hexa

/* CellSolution() for Cell3D_Hexa */
template<> inline
const double & ComputationalCell<ThreeD,Cell3D_Hexa,double>::
 CellSolution(int &) const {
  return U_cell;
}

template<> inline
double & ComputationalCell<ThreeD,Cell3D_Hexa,double>::
 CellSolution(int &){
  return U_cell;
}


/* SolutionAtCoordinates() for Cell3D_Hexa */
template<> inline
double ComputationalCell<ThreeD,Cell3D_Hexa,double>::
SolutionAtCoordinates(double & X_Coord, double & Y_Coord, double & Z_Coord, int VarPosition){
  return SolutionAtCoordinates(X_Coord, Y_Coord, Z_Coord);
}

/* CellMCC() for Cell3D_Hexa */
template<> inline
double & ComputationalCell<ThreeD,Cell3D_Hexa,double>::
 CellMCC(int &){
  return MCC;
}

template<> inline
const double & ComputationalCell<ThreeD,Cell3D_Hexa,double>::
 CellMCC(int &) const{
  return MCC;
}

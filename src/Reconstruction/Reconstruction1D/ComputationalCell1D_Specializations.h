// ComputationalCell1D_Specializations.h defines specializations of member functions of template ComputationalCell1D

#ifndef _COMPUTATIONALCELL1D_SPECIALIZATIONS_INCLUDED
#define _COMPUTATIONALCELL1D_SPECIALIZATIONS_INCLUDED

/*******************************************************
 * CLASS Template: ComputationalCell1D                 *
 * Specialization of member functions                  *
 ******************************************************/

// Specialization for SolutionType == double & GeometryType == Cell1D_Uniform/Cell1D_NonUniform

/* CellSolution() for Cell1D_NonUniform */
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

/* CellSolution() for Cell1D_Uniform */
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

/* CellMCC() for Cell1D_NonUniform */
template<> inline
double & ComputationalCell<OneD,Cell1D_NonUniform,double>::
 CellMCC(int &){
  return MCC;
}

template<> inline
const double & ComputationalCell<OneD,Cell1D_NonUniform,double>::
 CellMCC(int &) const{
  return MCC;
}

/* CellMCC() for Cell1D_Uniform */
template<> inline
double & ComputationalCell<OneD,Cell1D_Uniform,double>::
 CellMCC(int &){
  return MCC;
}

template<> inline
const double & ComputationalCell<OneD,Cell1D_Uniform,double>::
 CellMCC(int &) const{
  return MCC;
}


#endif

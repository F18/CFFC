// ComputationalCell2D_Implementations.h defines the implementation of member functions of template ComputationalCell2D

/*******************************************************
 * TEMPLATIZED CLASS: ComputationalCell2D              *
 * Implementation of member functions                  *
 ******************************************************/

// Constructor
template<class GeometryType,class SolutionType> inline
ComputationalCell<TwoD,GeometryType,SolutionType>::ComputationalCell():
  geom(),U_cell(0.0),TD(),TD_FirstOrder(),GeomCoeff(),SubgridSolution(),rings(0),RO(0),
     FinalOrder(0),ErrorRecL1Norm(0.0),ErrorRecL2Norm(0.0),ErrorRecMaxNorm(0.0), LimitedCell(OFF){
  NbSubgridPoints.reserve(TwoD);
}

// Copy Constructor
template<class GeometryType,class SolutionType> inline
ComputationalCell<TwoD,GeometryType,SolutionType>::ComputationalCell( const ComputationalCell & rhs):
  geom(rhs.geom), U_cell(rhs.U_cell),rings(rhs.rings),RO(rhs.RO),FinalOrder(rhs.FinalOrder),
     ErrorRecL1Norm(rhs.ErrorRecL1Norm), ErrorRecL2Norm(rhs.ErrorRecL2Norm), ErrorRecMaxNorm(rhs.ErrorRecMaxNorm),
     NbSubgridPoints(rhs.NbSubgridPoints),SubgridSolution(rhs.SubgridSolution), TD(rhs.TD){
}

// AdjustNumberOfSubgridPoints -> this function is called before the subgrid is generated
//                                but after the number of subgrid points is initialized
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::AdjustNumberOfSubgridPoints(void){

  /* Set the proper number of points in the x and y direction */
  /* X dir */
  if (NbSubgridPoints[0] < 3){
    NbSubgridPoints[0] = 3;
  } else if (NbSubgridPoints[0]%2 != 1){
    NbSubgridPoints[0] ++; //only odd numbers
  }
  /* Y dir */
  if (NbSubgridPoints[1] < 3){
    NbSubgridPoints[1] = 3;
  } else if (NbSubgridPoints[1]%2 != 1){
    NbSubgridPoints[1] ++; //only odd numbers
  }
}

// Generate cell in 2D;
template<class GeometryType,class SolutionType>
void ComputationalCell<TwoD,GeometryType,SolutionType>::GenerateCell(const int NbX, const int NbY, const int OrderOfRec){

  RO = OrderOfRec;
  FinalOrder = RO;
  NbSubgridPoints.reserve(TwoD);
  NbSubgridPoints.push_back(NbX);
  NbSubgridPoints.push_back(NbY);

  // adjust the number of subgrid points
  AdjustNumberOfSubgridPoints();

  // generate SubgridMesh
  SubgridSolution.SetMesh(NbSubgridPoints[0],NbSubgridPoints[1]);

  // generate the container of Taylor derivatives
  TD.GenerateContainer(OrderOfRec);
  TD_FirstOrder.GenerateContainer(1); /* OrderOfRec = 1 */

  // generate the container of the geometric coefficients
  GeomCoeff.GenerateContainer(OrderOfRec);
}

// Set Rings
template< class GeometryType, class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::SetRings(void){
  int ND;
  ND = TD.size();

  // TwoD
  switch(ND){
  case 1:  rings= 0; break;
  case 3:  rings= 1; break;
  case 6:  rings= 2; break;  // -> this is the default choice
  case 10: rings= 2; break;
  case 15: rings= 2; break;  //its different than the 1D situation (rings=3 would be the same as in 1D)
  case 21: rings= 3; break;
  case 28: rings= 3; break;
    /*       case 36: rings= 3; break; */
    /*       case 45: rings= 4; break; */
    /*       case 55: rings= 4; break; */
  default:
    rings = -1;
  }
}

// SetCell Parameters
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
SetCell (const GeometryType & geom_, const std::vector<int> & SubGridPoints, const int & RO_){
  require(TwoD <= SubGridPoints.size(), "ComputationalCell::SetCell() failed. The size of the provided vector is too small for the required dimension\n");
  // copy geometry
  SetCell(geom_);

  // generate SubgridMesh, TaylorDerivatives and Geometric Coefficients containers
  GenerateCell(SubGridPoints[0], SubGridPoints[1], RO_);

  // compute the geometric coefficients
  ComputeGeometricCoefficients();

  // Set number of rings
  SetRings();

  // Set geometry of the subgrid
  SetSubGridGeometry();
}

template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  SetCell (const std::vector<int> &SubGridPoints, const int & RO_){
  require(TwoD <= SubGridPoints.size(), "ComputationalCell::SetCell() failed. The size of the provided vector is too small for the required dimension\n");

  // generate SubgridMesh, TaylorDerivatives and Geometric Coefficients containers
  GenerateCell(SubGridPoints[0], SubGridPoints[1], RO_);

  // Set number of rings
  SetRings();
}

// SetSubGridGeometry()
/******************************************************************
* The subgrid is generated on a unit square and then is projected *
* in the physical space using a bilinear transformation.          *
*                                                                 *
* The coefficients of the transformation are obtained using the   *
* node of the cell.                                               *
*         x = a0 + a1*p + a2*q + a3*p*q                           *
*         y = b0 + b1*p + b2*q + b3*p*q                           *
******************************************************************/
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::SetSubGridGeometry(){

  double a0,a1,a2,a3;
  double b0,b1,b2,b3;

  /* Determine the coefficients of the bilinear transformation */
  a0 = geom.NodeSW().x;
  a1 = geom.NodeSE().x - geom.NodeSW().x;
  a2 = geom.NodeNW().x - geom.NodeSW().x;
  a3 = geom.NodeSW().x + geom.NodeNE().x - geom.NodeNW().x - geom.NodeSE().x;

  b0 = geom.NodeSW().y;
  b1 = geom.NodeSE().y - geom.NodeSW().y;
  b2 = geom.NodeNW().y - geom.NodeSW().y;
  b3 = geom.NodeSW().y + geom.NodeNE().y - geom.NodeNW().y - geom.NodeSE().y;

  double XdeltaSG, YdeltaSG;
  Vector nodePQ;			/* node in the computational plane */
  Vector nodeXY;			/* node in the physical plane */

  XdeltaSG = 1.0/(NbSubgridPoints[0]-1);
  YdeltaSG = 1.0/(NbSubgridPoints[1]-1);

  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j){
      /* Generate the points on the unit square */
      nodePQ.x = i*XdeltaSG;
      nodePQ.y = j*YdeltaSG;
      /* Transform in the physical plane */
      nodeXY.x = a0 + a1*nodePQ.x + a2*nodePQ.y + a3*nodePQ.x*nodePQ.y;
      nodeXY.y = b0 + b1*nodePQ.x + b2*nodePQ.y + b3*nodePQ.x*nodePQ.y;
      /* Write the node in the container */
      SubgridSolution(i,j).SetNode(nodeXY);
    }
}

// ComputeGeometricCoefficients()
/******************************************************************
 * Computes the geometric coefficients by computing integral      *
 * over the volume of the cell of the polynomial function having  *
 * the form                                                       *
 *       ((x-xCC)^n)*((y-yCC)^m) and divide by aria A             * 
******************************************************************/
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>
::ComputeGeometricCoefficients( ){

  int n,m;
  double DummyParam;

  /* Create the polynomial function for the cell */
  GeneralizedPolynomialFunctionOfTwoVariables Polynom(0,0,geom.xc.x,geom.xc.y);

#ifdef __Use_Iterator__
  /* Create iterator */
  GeomCoeffIterator Iter = GeomCoeff;

  for (Iter = GeomCoeff.begin(); Iter!=GeomCoeff.end(); Iter++){
    /* Compute the geometric coefficient */
    n = Iter->P1();
    m = Iter->P2();
    Polynom.ChangePowersTo(n,m);
    Iter->D() = IntegrateOverTheCell(Polynom,14,DummyParam)/geom.A();
  }
#else


  /* When the method of integration is used the moments associated with the first order are not quite accurate */
  for (int p1 = GeomCoeff.FirstElem(); p1<=GeomCoeff.LastElem(); ++p1){
    /* Compute the geometric coefficient */
    Polynom.ChangePowersTo(GeomCoeff(p1,true,true,true).P1(),GeomCoeff(p1,true,true,true).P2());
    GeomCoeff(p1,true,true,true).D() = IntegrateOverTheCell(Polynom,14,DummyParam)/geom.A();
  }

  /***************************************************************************** 
   * Subroutine for computing the geometric moments up to 3rd order            *
   * with respect to the cell centroid(xCC,yCC).                               *
   * These moments are integrals over the domain of the cell of the            *
   * polynomial function having the form ((x-xCC)^n)*((y-yCC)^m) and divided   *
   * by aria A.                                                                *
   ****************************************************************************/

  /* Moments for quadrilateral cells with straight boundaries */

  // Compute the coefficients of the bilinear interpolation for transforming
  // a quadrilatral cell into a unit lenght square

#if 0
  /* determine the coefficients of the transformation */
  double a0C = CellGeometry().NodeSW().x - CellGeometry().xc.x;
  double a1  = CellGeometry().NodeSE().x - CellGeometry().NodeSW().x;
  double a2  = CellGeometry().NodeNW().x - CellGeometry().NodeSW().x;
  double a3  = CellGeometry().NodeSW().x + CellGeometry().NodeNE().x - CellGeometry().NodeNW().x - CellGeometry().NodeSE().x;
  
  double b0C = CellGeometry().NodeSW().y - CellGeometry().xc.y;
  double b1  = CellGeometry().NodeSE().y - CellGeometry().NodeSW().y;
  double b2  = CellGeometry().NodeNW().y - CellGeometry().NodeSW().y;
  double b3  = CellGeometry().NodeSW().y + CellGeometry().NodeNE().y - CellGeometry().NodeNW().y - CellGeometry().NodeSE().y;

  double J0  = a1*b2 - a2*b1;
  double J1  = a1*b3 - a3*b1;
  double J2  = b2*a3 - a2*b3;


  switch(CellRecOrder()){

  case 3:   // Cubic moments
    // coefficient (3,0)
    GeomCoeff(3,0) = ((240*J0 + 120*J1 + 120*J2)*a0C*a0C*a0C + (360*J0 + 240*J1 + 180*J2)*a0C*a0C*a1  + 
		      (360*J0 + 180*J1 + 240*J2)*a0C*a0C*a2  + (180*J0 + 120*J1 + 120*J2)*a0C*a0C*a3  + 
		      (240*J0 + 180*J1 + 120*J2)*a0C*a1*a1  + (360*J0 + 240*J1 + 240*J2)*a0C*a1*a2  + 
		      (240*J0 + 180*J1 + 160*J2)*a0C*a1*a3  + (240*J0 + 120*J1 + 180*J2)*a0C*a2*a2  + 
		      (240*J0 + 160*J1 + 180*J2)*a0C*a2*a3  + (80*J0 + 60*J1 + 60*J2)*a0C*a3*a3  + 
		      (60*J0 + 48*J1 + 30*J2)*a1*a1*a1  + (120*J0 + 90*J1 + 80*J2)*a1*a1*a2  + 
		      (90*J0 + 72*J1 + 60*J2)*a1*a1*a3  + (120*J0 + 80*J1 + 90*J2)*a1*a2*a2  + 
		      (160*J0 + 120*J1 + 120*J2)*a1*a2*a3  + (60*J0 + 48*J1 + 45*J2)*a1*a3*a3  + 
		      (60*J0 + 30*J1 + 48*J2)*a2*a2*a2 + (90*J0 + 60*J1 + 72*J2)*a2*a2*a3  + 
		      (60*J0 + 45*J1 + 48*J2)*a2*a3*a3  + (15*J0 + 12*J1 + 12*J2)*a3*a3*a3) /geom.A() /0.240e3;

    // coefficient (2,1)
    GeomCoeff(2,1) = (((360*J0 + 180*J1 + 240*J2)*b2*a0C*a0C)  + ((360*J0 + 240*J1 + 240*J2)*b2*a0C*a1)  + 
		      ((480*J0 + 240*J1 + 360*J2)*b2*a0C*a2)  + ((240*J0 + 160*J1 + 180*J2)*b2*a0C*a3)  + 
		      ((120*J0 + 90*J1 + 80*J2)*b2*a1*a1)  + ((240*J0 + 160*J1 + 180*J2)*b2*a1*a2)  + 
		      ((160*J0 + 120*J1 + 120*J2)*b2*a1*a3)  + ((180*J0 + 90*J1 + 144*J2)*b2*a2*a2)  + 
		      ((180*J0 + 120*J1 + 144*J2)*b2*a2*a3)  + ((60*J0 + 45*J1 + 48*J2)*b2*a3*a3)  + 
		      ((720*J0 + 360*J1 + 360*J2)*a0C*a0C*b0C)  + ((360*J0 + 240*J1 + 180*J2)*a0C*a0C*b1)  + 
		      ((180*J0 + 120*J1 + 120*J2)*a0C*a0C*b3)  + ((720*J0 + 480*J1 + 360*J2)*a0C*a1*b0C)  + 
		      ((480*J0 + 360*J1 + 240*J2)*a0C*a1*b1)  + ((240*J0 + 180*J1 + 160*J2)*a0C*a1*b3)  + 
		      ((720*J0 + 360*J1 + 480*J2)*a0C*a2*b0C)  + ((360*J0 + 240*J1 + 240*J2)*a0C*a2*b1)  + 
		      ((240*J0 + 160*J1 + 180*J2)*a0C*a2*b3)  + ((360*J0 + 240*J1 + 240*J2)*a0C*a3*b0C)  + 
		      ((240*J0 + 180*J1 + 160*J2)*a0C*a3*b1)  + ((160*J0 + 120*J1 + 120*J2)*a0C*a3*b3)  + 
		      ((240*J0 + 180*J1 + 120*J2)*a1*a1*b0C)  + ((180*J0 + 144*J1 + 90*J2)*a1*a1*b1)  + 
		      ((90*J0 + 72*J1 + 60*J2)*a1*a1*b3)  + ((360*J0 + 240*J1 + 240*J2)*a1*a2*b0C)  + 
		      ((240*J0 + 180*J1 + 160*J2)*a1*a2*b1)  + ((160*J0 + 120*J1 + 120*J2)*a1*a2*b3)  + 
		      ((240*J0 + 180*J1 + 160*J2)*a1*a3*b0C)  + ((180*J0 + 144*J1 + 120*J2)*a1*a3*b1)  + 
		      ((120*J0 + 96*J1 + 90*J2)*a1*a3*b3)  + ((240*J0 + 120*J1 + 180*J2)*a2*a2*b0C)  + 
		      ((120*J0 + 80*J1 + 90*J2)*a2*a2*b1)  + ((90*J0 + 60*J1 + 72*J2)*a2*a2*b3)  + 
		      ((240*J0 + 160*J1 + 180*J2)*a2*a3*b0C)  + ((160*J0 + 120*J1 + 120*J2)*a2*a3*b1)  + 
		      ((120*J0 + 90*J1 + 96*J2)*a2*a3*b3)  + ((80*J0 + 60*J1 + 60*J2)*a3*a3*b0C)  + 
		      ((60*J0 + 48*J1 + 45*J2)*a3*a3*b1) + (45*J0 + 36*J1 + 36*J2)*a3*a3*b3) /geom.A() /0.720e3;
    
    // coefficient (1,2)
    GeomCoeff(1,2) = (((240*J0 + 120*J1 + 180*J2)*b2*b2*a0C) + ((120*J0 + 80*J1 + 90*J2)*b2*b2*a1) + 
		      ((180*J0 + 90*J1 + 144*J2)*b2*b2*a2) + ((90*J0 + 60*J1 + 72*J2)*b2*b2*a3) + 
		      ((720*J0 + 360*J1 + 480*J2)*b2*a0C*b0C) + ((360*J0 + 240*J1 + 240*J2)*b2*a0C*b1) + 
		      ((240*J0 + 160*J1 + 180*J2)*b2*a0C*b3) + ((360*J0 + 240*J1 + 240*J2)*b2*a1*b0C) + 
		      ((240*J0 + 180*J1 + 160*J2)*b2*a1*b1) + ((160*J0 + 120*J1 + 120*J2)*b2*a1*b3) + 
		      ((480*J0 + 240*J1 + 360*J2)*b2*a2*b0C) + ((240*J0 + 160*J1 + 180*J2)*b2*a2*b1) + 
		      ((180*J0 + 120*J1 + 144*J2)*b2*a2*b3) + ((240*J0 + 160*J1 + 180*J2)*b2*a3*b0C) + 
		      ((160*J0 + 120*J1 + 120*J2)*b2*a3*b1) + ((120*J0 + 90*J1 + 96*J2)*b2*a3*b3) + 
		      ((720*J0 + 360*J1 + 360*J2)*a0C*b0C*b0C) + ((720*J0 + 480*J1 + 360*J2)*a0C*b0C*b1) + 
		      ((360*J0 + 240*J1 + 240*J2)*a0C*b0C*b3) + ((240*J0 + 180*J1 + 120*J2)*a0C*b1*b1) + 
		      ((240*J0 + 180*J1 + 160*J2)*a0C*b1*b3) + ((80*J0 + 60*J1 + 60*J2)*a0C*b3*b3) + 
		      ((360*J0 + 240*J1 + 180*J2)*a1*b0C*b0C) + ((480*J0 + 360*J1 + 240*J2)*a1*b0C*b1) + 
		      ((240*J0 + 180*J1 + 160*J2)*a1*b0C*b3) + ((180*J0 + 144*J1 + 90*J2)*a1*b1*b1) + 
		      ((180*J0 + 144*J1 + 120*J2)*a1*b1*b3) + ((60*J0 + 48*J1 + 45*J2)*a1*b3*b3) + 
		      ((360*J0 + 180*J1 + 240*J2)*a2*b0C*b0C) + ((360*J0 + 240*J1 + 240*J2)*a2*b0C*b1) + 
		      ((240*J0 + 160*J1 + 180*J2)*a2*b0C*b3) + ((120*J0 + 90*J1 + 80*J2)*a2*b1*b1) + 
		      ((160*J0 + 120*J1 + 120*J2)*a2*b1*b3) + ((60*J0 + 45*J1 + 48*J2)*a2*b3*b3) + 
		      ((180*J0 + 120*J1 + 120*J2)*a3*b0C*b0C) + ((240*J0 + 180*J1 + 160*J2)*a3*b0C*b1) + 
		      ((160*J0 + 120*J1 + 120*J2)*a3*b0C*b3) + ((90*J0 + 72*J1 + 60*J2)*a3*b1*b1) + 
		      ((120*J0 + 96*J1 + 90*J2)*a3*b1*b3) + (45*J0 + 36*J1 + 36*J2)*a3*b3*b3) /geom.A() /0.720e3;

    // coefficient (0,3)
    GeomCoeff(0,3) = ((60*J0 + 30*J1 + 48*J2)*b2*b2*b2 + (240*J0 + 120*J1 + 180*J2)*b2*b2*b0C +
		      (120*J0 + 80*J1 + 90*J2)*b2*b2*b1 + (90*J0 + 60*J1 + 72*J2)*b2*b2*b3 + 
		      (360*J0 + 180*J1 + 240*J2)*b2*b0C*b0C + (360*J0 + 240*J1 + 240*J2)*b2*b0C*b1 + 
		      (240*J0 + 160*J1 + 180*J2)*b2*b0C*b3 + (120*J0 + 90*J1 + 80*J2)*b2*b1*b1 +
		      (160*J0 + 120*J1 + 120*J2)*b2*b1*b3 + (60*J0 + 45*J1 + 48*J2)*b2*b3*b3 + 
		      (240*J0 + 120*J1 + 120*J2)*b0C*b0C*b0C + (360*J0 + 240*J1 + 180*J2)*b0C*b0C*b1 + 
		      (180*J0 + 120*J1 + 120*J2)*b0C*b0C*b3 + (240*J0 + 180*J1 + 120*J2)*b0C*b1*b1 +
		      (240*J0 + 180*J1 + 160*J2)*b0C*b1*b3 + (80*J0 + 60*J1 + 60*J2)*b0C*b3*b3 + 
		      (60*J0 + 48*J1 + 30*J2)*b1*b1*b1 + (90*J0 + 72*J1 + 60*J2)*b1*b1*b3 + 
		      (60*J0 + 48*J1 + 45*J2)*b1*b3*b3 + (15*J0 + 12*J1 + 12*J2)*b3*b3*b3) /geom.A() /0.240e3;

  case 2:   // Quadratic moments
    // coefficient (2,0)
    GeomCoeff(2,0) = (((36*J0 + 18*J1 + 18*J2)*a0C*a0C) + ((36*J0 + 24*J1 + 18*J2)*a0C*a1) + 
		      ((36*J0 + 18*J1 + 24*J2)*a0C*a2) + ((18*J0 + 12*J1 + 12*J2)*a0C*a3) +
		      ((12*J0 + 9*J1 + 6*J2)*a1*a1) + ((18*J0 + 12*J1 + 12*J2)*a1*a2) + 
		      ((12*J0 + 9*J1 + 8*J2)*a1*a3) + ((12*J0 + 6*J1 + 9*J2)*a2*a2) + 
		      ((12*J0 + 8*J1 + 9*J2)*a2*a3) + ((4*J0 + 3*J1 + 3*J2)*a3*a3)) /geom.A() /0.36e2;
    
    // coefficient (0,2)
    GeomCoeff(0,2) = (((12*J0 + 6*J1 + 9*J2)*b2*b2) + ((36*J0 + 18*J1 + 24*J2)*b2*b0C) +
		      ((18*J0 + 12*J1 + 12*J2)*b2*b1) + ((12*J0 + 8*J1 + 9*J2)*b2*b3) + 
		      ((36*J0 + 18*J1 + 18*J2)*b0C*b0C) + ((36*J0 + 24*J1 + 18*J2)*b0C*b1) + 
		      ((18*J0 + 12*J1 + 12*J2)*b0C*b3) + ((12*J0 + 9*J1 + 6*J2)*b1*b1) + 
		      ((12*J0 + 9*J1 + 8*J2)*b1*b3) + ((4*J0 + 3*J1 + 3*J2)*b3*b3)) /geom.A() /0.36e2;
    
    // coefficient (1,1)
    GeomCoeff(1,1) = (((36*J0 + 18*J1 + 24*J2)*b2*a0C) + ((18*J0 + 12*J1 + 12*J2)*b2*a1) + 
		      ((24*J0 + 12*J1 + 18*J2)*b2*a2) + ((12*J0 + 8*J1 + 9*J2)*b2*a3) + 
		      ((72*J0 + 36*J1 + 36*J2)*a0C*b0C) + ((36*J0 + 24*J1 + 18*J2)*a0C*b1) +
		      ((18*J0 + 12*J1 + 12*J2)*a0C*b3) + ((36*J0 + 24*J1 + 18*J2)*a1*b0C) + 
		      ((24*J0 + 18*J1 + 12*J2)*a1*b1) + ((12*J0 + 9*J1 + 8*J2)*a1*b3) +
		      ((36*J0 + 18*J1 + 24*J2)*a2*b0C) + ((18*J0 + 12*J1 + 12*J2)*a2*b1) +
		      ((12*J0 + 8*J1 + 9*J2)*a2*b3) + ((18*J0 + 12*J1 + 12*J2)*a3*b0C) + 
		      ((12*J0 + 9*J1 + 8*J2)*a3*b1) + ((8*J0 + 6*J1 + 6*J2)*a3*b3)) / geom.A() /0.72e2;

  case 1:   // Linear moments
    GeomCoeff(1,0) = 0.0;
    GeomCoeff(0,1) = 0.0;
    
  case 0:   // Zero moment 
    GeomCoeff(0,0) = 1.0;
    
  }
#endif

#endif
};

// UpdateSubgridSolution()
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::UpdateSubgridSolution(){

  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j){
      SubgridSolution(i,j).SetValue( SolutionAt(SubgridSolution(i,j).GetNode()) );
    }
}

// template<typename FO>
// UpdateExactSubgridSolution(const FO)
template<class GeometryType,class SolutionType>
  template<typename FO> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  UpdateExactSubgridSolution(const FO FuncObj)
{
  for (int i=0; i<NbSubgridPoints[0]; ++i)
    for (int j=0; j<NbSubgridPoints[1]; ++j){
      SubgridSolution(i,j).SetValue( ExactSolutionAt(SubgridSolution(i,j).GetNode(), FuncObj) );
    }
}

// Compute solution : SolutionAt()
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<TwoD,GeometryType,SolutionType>::
  SolutionAt(double & deltaX, double & deltaY){
  return TD.ComputeSolutionFor(deltaX,deltaY);
}

//SolutionAtCoordinates(double & , double & )
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<TwoD,GeometryType,SolutionType>::
  SolutionAtCoordinates(double & X_Coord, double & Y_Coord){

  double DifferenceX, DifferenceY;
  DifferenceX = X_Coord - geom.xc.x;
  DifferenceY = Y_Coord - geom.xc.y;

  return TD.ComputeSolutionFor(DifferenceX,DifferenceY);
}

template<class GeometryType,class SolutionType> inline
double ComputationalCell<TwoD,GeometryType,SolutionType>::
  SolutionAtCoordinates(double & X_Coord, double & Y_Coord, int VarPosition){
  return SolutionAtCoordinates(X_Coord,Y_Coord)[VarPosition];
}

//SolutionAt(Node & )
template<class GeometryType,class SolutionType> inline
SolutionType ComputationalCell<TwoD,GeometryType,SolutionType>::
  SolutionAt(Node & node){

  Vector difference;

  difference = node.X - geom.xc;
  return TD.ComputeSolutionFor(difference.x,difference.y);
}

// CellSolution(int & VarPosition) const;
template<class GeometryType,class SolutionType> inline
const double & ComputationalCell<TwoD,GeometryType,SolutionType>::
  CellSolution(int & VarPosition) const {

  return U_cell[VarPosition];
}

// CellSolution(int & VarPosition);
template<class GeometryType,class SolutionType> inline
double & ComputationalCell<TwoD,GeometryType,SolutionType>::
  CellSolution(int & VarPosition) {
  
  return U_cell[VarPosition];
}

// CellMCC(int & VarPosition) const;
template<class GeometryType,class SolutionType> inline
const double & ComputationalCell<TwoD,GeometryType,SolutionType>::
  CellMCC(int & VarPosition) const {

  return MCC[VarPosition];
}

// CellMCC(int & VarPosition);
template<class GeometryType,class SolutionType> inline
double & ComputationalCell<TwoD,GeometryType,SolutionType>::
  CellMCC(int & VarPosition) {
  
  return MCC[VarPosition];
}

// template<typename PointerToTestFunction>
// SetInitialSolution(PointerToTestFunction)
/***************************************************************
Sets the value of the U_cell based on the integration of PtrFunc
over the domain of the cell
*****************************************************************/
template<class GeometryType,class SolutionType>
  template<typename PointerToTestFunction> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::SetInitialSolution(PointerToTestFunction PtrFunc)
{
  /* compute the average solution */
  U_cell = IntegrateOverTheCell(PtrFunc,6,U_cell)/geom.A();
  /* update the first coefficient of the Taylor expansion */
  TD(0,0) = U_cell;
  TD_FirstOrder(0,0) = U_cell;
}

// template<typename FunctionObject>
// void ComputeReconstructionError(const FO FuncObj)
/*******************************************************************
Computes the error of the reconstructed function in all three norms
********************************************************************/
template<class GeometryType,class SolutionType>
  template<typename FO>
void ComputationalCell<TwoD,GeometryType,SolutionType>::ComputeReconstructionError(const FO FuncObj){

  /* construct the error function for the L1 norm*/
  _Error_<FO,CompCellType, SolutionType> ErrorFunction(FuncObj,this);
  /* integrate the error */
  ErrorRecL1Norm = IntegrateOverTheCell(ErrorFunction,8,ErrorRecL1Norm);

  /* construct the error function for the L2 norm*/
  _Error_Square_<FO,CompCellType, SolutionType> ErrorFunctionSquare(FuncObj,this);
  /* integrate the error */
  ErrorRecL2Norm = IntegrateOverTheCell(ErrorFunctionSquare,8,ErrorRecL2Norm);

  /* get the maximum error based on the error values at the subgrid points */
  ErrorRecMaxNorm = SolutionType(0.0);
  for (int i=0; i<NbSubgridPoints[0]; ++i){
    for (int j=0; j<NbSubgridPoints[1]; ++j){
      ErrorRecMaxNorm = max(ErrorRecMaxNorm,ErrorFunction(SubgridSolution(i,j).GetNode().x(), 
							  SubgridSolution(i,j).GetNode().y()));
    }
  }
}

// OutputSubgridTecplot()
/* Outputs the coordinates of the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
   OutputMeshNodesTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = \"" << " Mesh nodes \"\n";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n";
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "F = POINT \n";
    output_file << "DT = (DOUBLE DOUBLE) \n";
  }

  output_file << CellGeometry().NodeSW() << "\n";

  output_file << setprecision(6);
}


// OutputSubgridTecplot()
/* Outputs the coordinates of the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputSubgridTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell2D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n"
		<< "ZONE T = \" Cell [" << iCell <<"," << jCell << "] \" \\ \n";
  } else {
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell << "] \" \\ \n";
  }

  SubgridSolution.OutputGeometryTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionTecplot()
/* Outputs the solution at the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputSolutionTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell2D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \n";
    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell << "] \" \\ \n";

  } else {
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell << "] \" \\ \n";
  }

  SubgridSolution.OutputSolutionTecplot(output_file);

  output_file << setprecision(6);
}

// OutputSolutionCellTecplot()
/* Outputs the solution at the centroid of the cell */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputSolutionCellTecplot(std::ofstream &output_file,const int & iCell,const int & jCell,const bool Title)
{
  output_file << setprecision(14);
  if (Title) {
    output_file << "TITLE = ComputationalCell2D (Subgrid Node Locations)\n ";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n";
    VarNames.PrintHeaderTecplot(output_file);
    output_file << "ZONE T = \" Cell [" << iCell <<"," << jCell << "] \" \\ \n";
    VarNames.PrintDataTypeTecplot(output_file,TwoD);
  }

  output_file << geom.xc << "\t" << TD(0,0) << "\n";

  output_file << setprecision(6);
}

// OutputSolutionCellTecplotOneZone()
/* Outputs the solution at the nodes of the subgrid for the specified SubgridNy */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const int &SubgridNy)
{
  SubgridSolution.OutputSolutionTecplotOneZone(output_file,SubgridNy);
}

// OutputPWC()
/* Outputs the piecewise constant solution at the nodes of the subgrid */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputPWC(std::ofstream &output_file, const int &SubgridNy)
{

  output_file << setprecision(14);

  for(int i=0; i<iSubgridPoints(); ++i)
    output_file << SubgridSolution(i,SubgridNy).GetNode() <<  "\t" <<  CellSolution() << "\n";

  output_file << setprecision(6);
}

// OutputSmoothnessIndicatorAtNodesTecplot()
/* Outputs the smoothness indicator at the nodes of the cell */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
  OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const int &SubgridNy)
{
  output_file << setprecision(14);

  for(int i=0; i<SubgridSolution.XPoints(); ++i)
    output_file << SubgridSolution(i,SubgridNy).GetNode() <<  "\t" << UnfitReconstructionFlag()  << "\n";

  output_file << setprecision(6);
}

// OutputSolutionAtGaussPointsTecplot()
/* Outputs the solution at the Gauss quadrature points of the cell */
template<class GeometryType,class SolutionType> inline
void ComputationalCell<TwoD,GeometryType,SolutionType>::
   OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file,  int Face, int GaussPoint)
{
  output_file << setprecision(14);

  /* Gauss Quadrature Point values */
  Vector2D GP1, GP2;

#if 0
  switch(Face){
  case SOUTH:
    CellGeometry().GaussQuadPointsFaceS(GP1,GP2);
    output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n"
		<< GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";
    break;

  case NORTH:
    CellGeometry().GaussQuadPointsFaceN(GP1,GP2);
    output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n"
		<< GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";
    break;

  case WEST:
    if (GaussPoint == 1){
      CellGeometry().GaussQuadPointsFaceW(GP1,GP2);
      output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";

      CellGeometry().GaussQuadPointsFaceE(GP1,GP2);
      output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";
    } else {
      CellGeometry().GaussQuadPointsFaceW(GP1,GP2);
      output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";

      CellGeometry().GaussQuadPointsFaceE(GP1,GP2);
      output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";
    }
    break;
  }
#endif


  switch(Face){
  case 1:
    // W2
    CellGeometry().GaussQuadPointsFaceW(GP1,GP2);
    output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";

    // S1 & S2
    CellGeometry().GaussQuadPointsFaceS(GP1,GP2);
    output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n"
		<< GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";

    // E1
    CellGeometry().GaussQuadPointsFaceE(GP1,GP2);
    output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";
    break;

  case 2:
    // W1
    CellGeometry().GaussQuadPointsFaceW(GP1,GP2);
    output_file << GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";

    // N2 & N1
    CellGeometry().GaussQuadPointsFaceN(GP1,GP2);
    output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n"
		<< GP1 << "\t" << SolutionAtCoordinates(GP1.x,GP1.y) << "\n";

    // E2
    CellGeometry().GaussQuadPointsFaceE(GP1,GP2);
    output_file << GP2 << "\t" << SolutionAtCoordinates(GP2.x,GP2.y) << "\n";
    break;
  }

  output_file << setprecision(6);
}

// Friend functions
// Operator ==
template<class GeometryType,class SolutionType> inline
bool operator==(const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
		const ComputationalCell<TwoD,GeometryType,SolutionType>& right){

  return (left.geom==right.geom)&&(left.U_cell==right.U_cell)&&(left.TD==right.TD)&&(left.NbSubgridPoints==right.NbSubgridPoints);
}

// Operator !=
template<class GeometryType,class SolutionType> inline
bool operator!=(const ComputationalCell<TwoD,GeometryType,SolutionType>& left,
		const ComputationalCell<TwoD,GeometryType,SolutionType>& right){
  return !(left==right);
}

// Operator <<
template<class GeometryType,class SolutionType> inline
std::ostream& operator<< (std::ostream& os, const ComputationalCell<TwoD,GeometryType,SolutionType>& Obj){
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
HeaderData ComputationalCell<TwoD,GeometryType,SolutionType>::VarNames("Solution");

/* Include full specialization of template member functions */
#include "ComputationalCell2D_Specializations.h"

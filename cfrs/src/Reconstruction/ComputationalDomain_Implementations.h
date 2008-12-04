/* ComputationalDomain_Implementations.h: Header file defining the implementations of 
                                          template computational domain member functions */

#include "../../../src_2D/HighOrderReconstruction/ReconstructionHelpers.h"

/****************************************************************************
 * TEMPLATIZED CLASS: ComputationalDomain                                   *
 * Implementation of member functions.                                      *
 ***************************************************************************/

// Default Constructor
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ComputationalDomain()
  :SolnPtr(NULL), L1Norm(0.0), L2Norm(0.0), LMaxNorm(0.0)
{
  N_XYZ.reserve(SpaceDimension);
  IndexLow.reserve(SpaceDimension);
  IndexUp.reserve(SpaceDimension);
  Nghost = 0;
  for (int i=0; i<SpaceDimension; ++i){
    N_XYZ.push_back(0);
    IndexLow.push_back(0);
    IndexUp.push_back(0);
  }
  MaxDeltaSolutionOverDomain.reserve(NumberOfParameters);
  CharacteristicLength = ONE;
}

//Copy constructor
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ComputationalDomain(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType> & rhs){


  // copy the dimensions
  N_XYZ = rhs.N_XYZ;
  IndexLow = rhs.IndexLow;
  IndexUp = rhs.IndexUp;
  Nghost = rhs.Nghost;
  L1Norm = rhs.L1Norm;
  L2Norm = rhs.L2Norm;
  LMaxNorm = rhs.LMaxNorm;
  MaxDeltaSolutionOverDomain = rhs.MaxDeltaSolutionOverDomain;
  CharacteristicLength = rhs.CharacteristicLength;

  // allocate memory
  allocate();

  // copy the elements
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<N_XYZ[0]; ++i){
      SolnPtr[0][0][i] = rhs.SolnPtr[0][0][i];
    }
    break;
  case TwoD:
    for (int i=0; i<N_XYZ[0]; ++i)
      for (int j=0; j<N_XYZ[1]; ++j)
	SolnPtr[0][j][i] = rhs.SolnPtr[0][j][i];
    break;
  case ThreeD:
    for (int i=0; i<N_XYZ[0]; ++i)
      for (int j=0; j<N_XYZ[1]; ++j)
	for (int k=0; k<N_XYZ[2]; ++k)
	  SolnPtr[k][j][i] = rhs.SolnPtr[k][j][i];
    break;
  }
}

// Assignment operator
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  ComputationalDomain<SpaceDimension,GeometryType,SolutionType>  & 
  ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::operator=(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType> & rhs) { 
  
  if (this == &rhs) return *this;

  // if the object has memory allocated this is going to be deleted
  deallocate();

  // copy the dimensions
  N_XYZ = rhs.N_XYZ;
  IndexLow = rhs.IndexLow;
  IndexUp = rhs.IndexUp;
  Nghost = rhs.Nghost;
  L1Norm = rhs.L1Norm;
  L2Norm = rhs.L2Norm;
  LMaxNorm = rhs.LMaxNorm;
  MaxDeltaSolutionOverDomain = rhs.MaxDeltaSolutionOverDomain;
  CharacteristicLength = rhs.CharacteristicLength;

  // allocate memory
  allocate();
  // copy the elements
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<N_XYZ[0]; ++i)
      SolnPtr[0][0][i] = rhs.SolnPtr[0][0][i];
    break;
  case TwoD:
    for (int i=0; i<N_XYZ[0]; ++i)
      for (int j=0; j<N_XYZ[1]; ++j)
	SolnPtr[0][j][i] = rhs.SolnPtr[0][j][i];
    break;
  case ThreeD:
    for (int i=0; i<N_XYZ[0]; ++i)
      for (int j=0; j<N_XYZ[1]; ++j)
	for (int k=0; k<N_XYZ[2]; ++k)
	  SolnPtr[k][j][i] = rhs.SolnPtr[k][j][i];
    break;
  }
  return *this;
}

/* Memory management functions */
// Allocate function
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::allocate( ) {
  switch(SpaceDimension){
  case OneD:
    SolnPtr = new CompCellType** ;
    SolnPtr[0] = new CompCellType*;
    SolnPtr[0][0] = new CompCellType[N_XYZ[0]]( );
    break;
  case TwoD:
    SolnPtr = new CompCellType** ;
    SolnPtr[0] = new CompCellType* [N_XYZ[1]];
    for (int i=0; i<N_XYZ[1]; ++i)
      SolnPtr[0][i] = new CompCellType[N_XYZ[0]]( );
    break;
  case ThreeD:
    SolnPtr = new CompCellType** [N_XYZ[2]];
    for (int i=0; i<N_XYZ[2]; ++i)
      SolnPtr[i] = new CompCellType* [N_XYZ[1]];
    for (int i=0; i<N_XYZ[2]; ++i)
      for (int j=0; j<N_XYZ[1]; ++j)
	SolnPtr[i][j] = new CompCellType[N_XYZ[0]]( );
    break;
  default: 
    std::cout << "In ComputationalDomain<>::allocate() the number of dimensions are not properly defined.\n";
  }
}

// Deallocate function
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::deallocate( ) {
  if (SolnPtr!=NULL){
    switch(SpaceDimension){
    case OneD:
      delete [] SolnPtr[0][0];
      delete SolnPtr[0];
      delete SolnPtr;
      break;
    case TwoD:
      for (int i=0; i<N_XYZ[1]; ++i)
	delete [] SolnPtr[0][i];
      delete [] SolnPtr[0];
      delete SolnPtr;
      break;
    case ThreeD:
      for (int i=0; i<N_XYZ[2]; ++i)
	for (int j=0; j<N_XYZ[1]; ++j)
	  delete [] SolnPtr[i][j];
      for (int i=0; i<N_XYZ[2]; ++i)
	delete [] SolnPtr[i];
      delete [] SolnPtr;
      break;
    default: 
      std::cout << "In ComputationalDomain::deallocate() the number of dimensions are not properly defined.\n";
    }

    SolnPtr = NULL;
  }
}

template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
bool ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::null(void) const
{
  if (SolnPtr==NULL){
    return true;
  } else {
    return false;
  }
}

// SetDomain(int,int)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const int & NCx, const int & _Nghost_)
{
  require(N_XYZ.size()==OneD, "ComputationalDomain::SetDomain failed\n");
  deallocate();
  Nghost = _Nghost_;
  N_XYZ[0] = NCx + 2*Nghost;

  /* set indexes */
  IndexLow[0] = Nghost;
  IndexUp[0] = NCx + Nghost - 1;

  allocate();
}

// SetDomain(Reconstruct1D_Input_Parameters & IP) -> Set a 1D uniform mesh
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const Reconstruct1D_Input_Parameters & IP)
{

  std::cout << "    Allocate memory for the Computational Domain:\n";
  // allocate memory
  SetDomain(IP.iCell(), IP.Nghost());


  // set the neccessary local parameters
  vector<int> SubGridPoints;
  SubGridPoints.reserve(1);
  SubGridPoints.push_back( IP.NumSubGridI() );

  // Local variables
  double FirstCentroid, DeltaX;
  Cell1D_NonUniform InitCell;

  std::cout << "     Set the Geometry:\n";

  // Set the geometry
  switch(IP.i_Grid){
  case GRID_READ_FROM_GRID_DATA_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Data_File subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    break;
  case GRID_READ_FROM_DEFINITION_FILE:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Read_From_Definition_File subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    break;
  case GRID_NONUNIFORM:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cout << "\n Grid_Nonuniform subrutine wasn't implemented yet\a" << endl;
    cout.flush();
    break;
  default:
    // uniform grid
    DeltaX = IP.LengthX()*IP.ScaleX()/IP.iCell();
    FirstCentroid = IP.MinX()*IP.ScaleX() + 0.5*DeltaX + IP.ShiftX(); // the centroid of the first computational cell

    for (int i=iStart(), Counter=0; i<=iEnd(); ++i, ++Counter){
      InitCell.x = FirstCentroid + Counter*DeltaX + IP.ShiftX();
      InitCell.dx = DeltaX;
      SolnPtr[0][0][i].SetCell(InitCell,SubGridPoints,IP.RecOrder());
    }

    // Set the ghost cells
    for (int i=0; i<iStart(); ++i){
      // left ghost cells
      InitCell.x  = 2*IP.MinX() - SolnPtr[0][0][iStart()+IP.Nghost()-1-i].CellCenter();
      InitCell.dx = SolnPtr[0][0][iStart()+IP.Nghost()-1-i].CellDelta();
      SolnPtr[0][0][i].SetCell(InitCell,SubGridPoints,IP.RecOrder());

      // right ghost cells
      InitCell.x  = 2*IP.MaxX() - SolnPtr[0][0][iEnd()-i].CellCenter();
      InitCell.dx = SolnPtr[0][0][iEnd()-i].CellDelta();
      SolnPtr[0][0][iEnd()+1+i].SetCell(InitCell,SubGridPoints,IP.RecOrder());
    }

  }

  // Set the characteristic length of the geometry
  CharacteristicLength = IP.CharacteristicLength;
  CutoffKnob = IP.CutoffKnob();

  // Reset the L1, L2 and LMax norms
  L1Norm = SolutionType(0.0);
  L2Norm = SolutionType(0.0);
  LMaxNorm = SolutionType(0.0);

}

// SetDomain(int,int,int)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const int & NCx, const int & NCy, const int & _Nghost_)
{
   require(N_XYZ.size()==TwoD, "ComputationalDomain::SetDomain failed\n");
   deallocate();
   Nghost = _Nghost_;
   N_XYZ[0] = NCx + 2*Nghost;
   N_XYZ[1] = NCy + 2*Nghost;

  /* set indexes */
  IndexLow[0] = Nghost;
  IndexLow[1] = Nghost;
  IndexUp[0] = NCx + Nghost - 1;
  IndexUp[1] = NCy + Nghost - 1;

  allocate();
}

// SetDomain(int,int,int,int)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const int & NCx, const int & NCy, const int & NCz, const int & _Nghost_)
{
   require(N_XYZ.size()==ThreeD, "ComputationalDomain::SetDomain failed\n");
   deallocate();
   Nghost = _Nghost_;
   N_XYZ[0] = NCx + 2*Nghost;
   N_XYZ[1] = NCy + 2*Nghost;
   N_XYZ[2] = NCz + 2*Nghost;

  /* set indexes */
  IndexLow[0] = Nghost;
  IndexLow[1] = Nghost;
  IndexLow[2] = Nghost;
  IndexUp[0] = NCx + Nghost - 1;
  IndexUp[1] = NCy + Nghost - 1;
  IndexUp[2] = NCz + Nghost - 1;

  allocate();
}

// SetDomain(const Grid2D_Quad_Block, const Reconstruct2D_Input_Parameters)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(Grid2D_Quad_Block & Grid, const Reconstruct2D_Input_Parameters & IP)
{

  int i,j;

  // allocate memory
  SetDomain(IP.iCell(), IP.jCell(), IP.Nghost());

  // allocate all the neccessary memory for the computational cell
  vector<int> SubGridPoints;
  SubGridPoints.reserve(2);
  SubGridPoints.push_back( IP.NumSubGridI() );
  SubGridPoints.push_back( IP.NumSubGridJ() );

  // Set the geometry
  Cell2D ijCell;
  Cell2D_Quad InitCell;

  for (j = 0; j< Grid.NCj; ++j)
    for (i = 0; i< Grid.NCi; ++i){
      ijCell.setindex(i,j);
      InitCell.setnodes(Grid.nodeSW(ijCell).X.x,Grid.nodeSW(ijCell).X.y,
			Grid.nodeSE(ijCell).X.x,Grid.nodeSE(ijCell).X.y,
			Grid.nodeNE(ijCell).X.x,Grid.nodeNE(ijCell).X.y,
			Grid.nodeNW(ijCell).X.x,Grid.nodeNW(ijCell).X.y);
      SolnPtr[0][j][i].SetCell(InitCell,SubGridPoints,IP.RecOrder());  
    }

  // Set the characteristic length of the geometry
  CharacteristicLength = IP.CharacteristicLength;
  CutoffKnob = IP.CutoffKnob();

  // Reset the L1, L2 and LMax norms
  L1Norm = SolutionType(0.0);
  L2Norm = SolutionType(0.0);
  LMaxNorm = SolutionType(0.0);

}

// SetDomain(const Reconstruct2D_Input_Parameters)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const Reconstruct2D_Input_Parameters & IP)
{

  ifstream In_File;
  int i,j, Nghost, RecOrder;

  // Open data file.
  In_File.open(IP.Input_File_Name, ios::in);
  if (In_File.bad()) {
    std::cout << "SetDomain()::ERROR!!! Can't open the input data file " << IP.Input_File_Name << "\n";
    std::cout.flush();
    exit(1);
  }

  // Read the number of cells in the Xdir and in the Ydir. Read number of ghost cells and the reconstruction order
  In_File >> i >> j >> Nghost >> RecOrder;

  // allocate memory
  SetDomain(i, j, Nghost);

  // allocate all the neccessary memory for the computational cell
  vector<int> SubGridPoints;
  SubGridPoints.reserve(2);
  SubGridPoints.push_back( IP.NumSubGridI() );
  SubGridPoints.push_back( IP.NumSubGridJ() );

  double XnodeSW, YnodeSW, XnodeSE, YnodeSE;
  double XnodeNE, YnodeNE, XnodeNW, YnodeNW;
  double CentroidX, CentroidY;
  GeometricIntegrals GeomCoeff; GeomCoeff.GenerateContainer(RecOrder);

  // Set the geometry
  Cell2D_Quad InitCell;

  for (i = 0; i<= iLastCell(); ++i){
    for (j = 0; j<= jLastCell(); ++j){

      In_File >> XnodeSW >> YnodeSW >> XnodeSE >> YnodeSE >> XnodeNE >> YnodeNE >> XnodeNW >> YnodeNW;

      // Set the nodes
      InitCell.setnodes(XnodeSW,YnodeSW,
			XnodeSE,YnodeSE,
			XnodeNE,YnodeNE,
			XnodeNW,YnodeNW);
      SolnPtr[0][j][i].SetCell(InitCell,SubGridPoints,RecOrder);  

      In_File >> CentroidX >> CentroidY;

      // Set the centroid
      SolnPtr[0][j][i].CellGeometry().xc.x = CentroidX;
      SolnPtr[0][j][i].CellGeometry().xc.y = CentroidY;

      In_File >> GeomCoeff;

      // Copy geometric coefficients
      SolnPtr[0][j][i].CellGeomCoeff() = GeomCoeff;

      In_File.setf(ios::skipws);
      In_File >> SolnPtr[0][j][i].CellSolution();


      SolnPtr[0][j][i].CellDeriv(0,true,true,true) = SolnPtr[0][j][i].CellSolution();
      SolnPtr[0][j][i].UpdateSubgridSolution();
    }
  }

  In_File.close();

  // Reset the L1, L2 and LMax norms
  L1Norm = SolutionType(0.0);
  L2Norm = SolutionType(0.0);
  LMaxNorm = SolutionType(0.0);
}

// SetDomain(const Grid3D_Hexa_Block, const Reconstruct3D_Input_Parameters)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(Grid3D_Hexa_Block & Grid, const Reconstruct3D_Input_Parameters & IP)
{

  int i,j,k;

  // allocate memory
  SetDomain(IP.iCell(), IP.jCell(), IP.kCell(), IP.Nghost());

  // allocate all the neccessary memory for the computational cell
  vector<int> SubGridPoints;
  SubGridPoints.reserve(3);
  SubGridPoints.push_back( IP.NumSubGridI() );
  SubGridPoints.push_back( IP.NumSubGridJ() );
  SubGridPoints.push_back( IP.NumSubGridK() );

  // Set the geometry
  Cell3D ijkCell;
  Cell3D_Hexa InitCell;

  for(k=0; k<=kLastCell(); ++k){
    for(j=0; j<=jLastCell(); ++j){
      for(i=0; i<=iLastCell(); ++i){

	ijkCell.setindex(i,j,k);
	InitCell.setnodes(Grid.nodeSWTop(ijkCell).X.x,Grid.nodeSWTop(ijkCell).X.y,Grid.nodeSWTop(ijkCell).X.z,
			  Grid.nodeSETop(ijkCell).X.x,Grid.nodeSETop(ijkCell).X.y,Grid.nodeSETop(ijkCell).X.z,
			  Grid.nodeNETop(ijkCell).X.x,Grid.nodeNETop(ijkCell).X.y,Grid.nodeNETop(ijkCell).X.z,
			  Grid.nodeNWTop(ijkCell).X.x,Grid.nodeNWTop(ijkCell).X.y,Grid.nodeNWTop(ijkCell).X.z,
			  Grid.nodeSWBot(ijkCell).X.x,Grid.nodeSWBot(ijkCell).X.y,Grid.nodeSWBot(ijkCell).X.z,
			  Grid.nodeSEBot(ijkCell).X.x,Grid.nodeSEBot(ijkCell).X.y,Grid.nodeSEBot(ijkCell).X.z,
			  Grid.nodeNEBot(ijkCell).X.x,Grid.nodeNEBot(ijkCell).X.y,Grid.nodeNEBot(ijkCell).X.z,
			  Grid.nodeNWBot(ijkCell).X.x,Grid.nodeNWBot(ijkCell).X.y,Grid.nodeNWBot(ijkCell).X.z);
	SolnPtr[k][j][i].SetCell(InitCell,SubGridPoints,IP.RecOrder());
      
      } /* end for */
    } /* end for */
  } /* end for */

  // Set the characteristic length of the geometry
  CharacteristicLength = IP.CharacteristicLength;
  CutoffKnob = IP.CutoffKnob();

  // Reset the L1, L2 and LMax norms
  L1Norm = SolutionType(0.0);
  L2Norm = SolutionType(0.0);
  LMaxNorm = SolutionType(0.0);
 
}

// SetDomain(const Reconstruct3D_Input_Parameters)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetDomain(const Reconstruct3D_Input_Parameters & IP)
{
 
  ifstream In_File;
  int i,j,k, Nghost, RecOrder;

  // Open data file.
  In_File.open(IP.Input_File_Name, ios::in);
  if (In_File.bad()) {
    std::cout << "SetDomain()::ERROR!!! Can't open the input data file " << IP.Input_File_Name << "\n";
    std::cout.flush();
    exit(1);
  }

  // Read the number of cells in the Xdir and in the Ydir. Read number of ghost cells and the reconstruction order
  In_File >> i >> j >> k >> Nghost >> RecOrder;

  // allocate memory
  SetDomain(i, j, k, Nghost);

  // allocate all the neccessary memory for the computational cell
  vector<int> SubGridPoints;
  SubGridPoints.reserve(3);
  SubGridPoints.push_back( IP.NumSubGridI() );
  SubGridPoints.push_back( IP.NumSubGridJ() );
  SubGridPoints.push_back( IP.NumSubGridK() );

  double XnodeSWTop, YnodeSWTop, ZnodeSWTop, XnodeSETop, YnodeSETop, ZnodeSETop;
  double XnodeNETop, YnodeNETop, ZnodeNETop, XnodeNWTop, YnodeNWTop, ZnodeNWTop;
  double XnodeSWBot, YnodeSWBot, ZnodeSWBot, XnodeSEBot, YnodeSEBot, ZnodeSEBot;
  double XnodeNEBot, YnodeNEBot, ZnodeNEBot, XnodeNWBot, YnodeNWBot, ZnodeNWBot;
  double CentroidX, CentroidY, CentroidZ;
  GeometricIntegrals GeomCoeff; 

  GeomCoeff.GenerateContainer(RecOrder);

  // Set the geometry
  Cell3D_Hexa InitCell;

  for (i = 0; i<= iLastCell(); ++i){
    for (j = 0; j<= jLastCell(); ++j){
      for(k = 0; k<= kLastCell(); ++k){

	In_File >> XnodeSWTop >> YnodeSWTop >> ZnodeSWTop >> XnodeSETop >> YnodeSETop >> ZnodeSETop 
		>> XnodeNETop >> YnodeNETop >> ZnodeNETop >> XnodeNWTop >> YnodeNWTop >> ZnodeNWTop;
	In_File >> XnodeSWBot >> YnodeSWBot >> ZnodeSWBot >> XnodeSEBot >> YnodeSEBot >> ZnodeSEBot 
		>> XnodeNEBot >> YnodeNEBot >> ZnodeNEBot >> XnodeNWBot >> YnodeNWBot >> ZnodeNWBot;

	// Set the nodes
	InitCell.setnodes(XnodeSWTop,YnodeSWTop,ZnodeSWTop,
			  XnodeSETop,YnodeSETop,ZnodeSETop, 
			  XnodeNETop,YnodeNETop,ZnodeNETop,
			  XnodeNWTop,YnodeNWTop,ZnodeNWTop,
			  XnodeSWBot,YnodeSWBot,ZnodeSWBot,
			  XnodeSEBot,YnodeSEBot,ZnodeSEBot, 
			  XnodeNEBot,YnodeNEBot,ZnodeNEBot,
			  XnodeNWBot,YnodeNWBot,ZnodeNWBot);

	SolnPtr[k][j][i].SetCell(InitCell,SubGridPoints,RecOrder);  
	
	In_File >> CentroidX >> CentroidY >> CentroidZ;
	
	// Set the centroid
	SolnPtr[k][j][i].CellGeometry().xc.x = CentroidX;
	SolnPtr[k][j][i].CellGeometry().xc.y = CentroidY;
	SolnPtr[k][j][i].CellGeometry().xc.z = CentroidZ;
	
	In_File >> GeomCoeff;
	
	// Copy geometric coefficients
	SolnPtr[k][j][i].CellGeomCoeff() = GeomCoeff;
	
	In_File.setf(ios::skipws);
	In_File >> SolnPtr[k][j][i].CellSolution();
	
	
	SolnPtr[k][j][i].CellDeriv(0) = SolnPtr[k][j][i].CellSolution();
	SolnPtr[k][j][i].UpdateSubgridSolution();

      } /* end for */
    } /* end for */
  } /* end for */

  In_File.close();

  // Reset the L1, L2 and LMax norms
  L1Norm = SolutionType(0.0);
  L2Norm = SolutionType(0.0);
  LMaxNorm = SolutionType(0.0);
}

/* Overloaded operators */
// bracket ( int )
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
typename ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::CompCellType & 
  ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::operator( )(const int &i){
  require(SpaceDimension == OneD, "ComputationalDomain::operator( ) failed." 
	  " You are not working in 1D. Wrong number of parameters passed!.\n");
  require((i<N_XYZ[0])&&(i>=0),"ComputationalDomain::operator( ) failed." 
	  " The indices are out of bounds! \n");

  return SolnPtr[0][0][i];
}

// bracket ( int, int )
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
typename ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::CompCellType & 
  ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::operator( )(const int &i,const int &j){
  require(SpaceDimension == TwoD, "ComputationalDomain::operator( ) failed."
	  " You are not working in 2D. Wrong number of parameters passed!.\n");
  require((i<N_XYZ[0])&&(i>=0)&&(j<N_XYZ[1])&&(j>=0), "ComputationalDomain::operator( ) failed."
	  " The indices are out of bounds! \n");

  return SolnPtr[0][j][i];
}

// bracket ( int, int , int)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
typename ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::CompCellType & 
  ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::operator( )(const int &i,const int &j,const int &k){
  require(SpaceDimension == ThreeD, "ComputationalDomain::operator( ) failed."
	  " You are not working in 3D. Wrong number of parameters passed!.\n");

  return SolnPtr[k][j][i];
}

// Operation functions
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
int ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::NumberOfNodes(void) const{
  switch(SpaceDimension){
  case OneD:
    return (N_XYZ[0]-2*Nghost+1);
    break;
  case TwoD:
    return (N_XYZ[0]-2*Nghost+1)*(N_XYZ[1]-2*Nghost+1);
    break;
  case ThreeD:
    return (N_XYZ[0]-2*Nghost+1)*(N_XYZ[1]-2*Nghost+1)*(N_XYZ[2]-2*Nghost+1);
    break;
  }
  return 0;
}

// OutputNodesTecplot
/* Outputs the coordinates of the nodes including the ghost cells */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputNodesTecplot(std::ofstream &output_file, const bool Title) const
{

  int i,j,k, _dummy_parameter_(0); 
  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(i=0; i<=iLastCell(); ++i){
      SolnPtr[0][0][i].OutputSubgridTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }

    break;

  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    for(j=0; j<=jLastCell(); ++j)
      for(i=0; i<=iLastCell(); ++i){
	SolnPtr[0][j][i].OutputSubgridTecplot(output_file,i,j,_dummy_parameter_);
      }
    
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    for(k=0; k<=kLastCell(); ++k)
      for(j=0; j<=jLastCell(); ++j)
	for(i=0; i<=iLastCell(); ++i){
	  SolnPtr[k][j][i].OutputSubgridTecplot(output_file,i,j,k);
	}

    break; 
  }
}

// OutputMeshNodesTecplot
/* Outputs the coordinates of the interior mesh nodes */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputMeshNodesTecplot(std::ofstream &output_file,const bool Title ) const
{
  output_file << setprecision(14);
  int i,j,k, _dummy_parameter_(0);

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << iEnd() - iStart() + 2 << "\n"
		<< "F = POINT \n";
    
    for(i=iStart(); i<=iEnd()+1; ++i){
      SolnPtr[0][0][i].OutputMeshNodesTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << iEnd() - iStart() + 2 << "\n"
		<< "J=" << jEnd() - jStart() + 2 << "\n"
		<< "F = POINT \n";

    for(j=jStart(); j<=jEnd()+1; ++j)
      for(i=iStart(); i<=iEnd()+1; ++i){
	SolnPtr[0][j][i].OutputMeshNodesTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << iEnd() - iStart() + 2 << "\n"
		<< "J=" << jEnd() - jStart() + 2 << "\n"
		<< "K=" << kEnd() - kStart() + 2 << "\n"
		<< "F = POINT \n";

    for(k=kStart(); k<=kEnd()+1; ++k)
      for(j=jStart(); j<=jEnd()+1; ++j)
	for(i=iStart(); i<=iEnd()+1; ++i){
	  SolnPtr[k][j][i].OutputMeshNodesTecplot(output_file,i,j,k);
	}

    break; 
  }
}

// OutputMeshNodesTecplot
/* Outputs the coordinates of the interior mesh nodes */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputMeshCellsTecplot(std::ofstream &output_file,const bool Title ) const
{
  output_file << setprecision(14);
  int i,j,k, _dummy_parameter_(0);

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << N_XYZ[0] << "\n"
		<< "F = POINT \n";
    
    for(i=iStart(); i<=iEnd()+1; ++i){
      SolnPtr[0][0][i].OutputMeshCellsTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << N_XYZ[0] << "\n"
		<< "J=" << N_XYZ[1] << "\n"
		<< "F = POINT \n";

    for(j=jStart(); j<=jEnd()+1; ++j)
      for(i=iStart(); i<=iEnd()+1; ++i){
	SolnPtr[0][j][i].OutputMeshCellsTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    output_file << "ZONE T = \" Mesh \" \\ \n"
		<< "I=" << N_XYZ[0] << "\n"
		<< "J=" << N_XYZ[1] << "\n"
		<< "K=" << N_XYZ[2] << "\n"
		<< "F = POINT \n";

    for(k=kStart()-Nghost; k<=kEnd()+Nghost; ++k)
      for(j=jStart()-Nghost; j<=jEnd()+Nghost; ++j)
	for(i=iStart()-Nghost; i<=iEnd()+Nghost; ++i){
	  SolnPtr[k][j][i].OutputMeshCellsTecplot(output_file,i,j,k);
	}

    break; 
  }
}


// OutputSolutionNodesTecplot
/* Outputs the solution at the nodes without the ghost cells */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputSolutionNodesTecplot(std::ofstream &output_file, const bool Title) const
{

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(int i=iStart(), _dummy_parameter_=0; i<=iEnd(); ++i){
	SolnPtr[0][0][i].OutputSolutionTecplot(output_file,i,_dummy_parameter_, _dummy_parameter_);
    }
    break;
    
  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    for(int j=jStart(), _dummy_parameter_=0; j<=jEnd(); ++j)
      for(int i=iStart(); i<=iEnd(); ++i){
	SolnPtr[0][j][i].OutputSolutionTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    for(int k=kStart(), _dummy_parameter_=0; k<=kEnd(); ++k)
      for(int j=jStart(); j<=jEnd(); ++j)
	for(int i=iStart(); i<=iEnd(); ++i){
	  SolnPtr[k][j][i].OutputSolutionTecplot(output_file,i,j,_dummy_parameter_);
	}
    break;
  }
}

// OutputFullSolutionNodesTecplot
/* Outputs the solution at the nodes including the ghost cells */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputFullSolutionNodesTecplot(std::ofstream &output_file, const bool Title) const
{

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(int i=0, _dummy_parameter_=0; i<=iLastCell(); ++i){
      SolnPtr[0][0][i].OutputSolutionTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(int j=0, _dummy_parameter_=0; j<=jLastCell(); ++j)
      for(int i=0; i<=iLastCell(); ++i){
	SolnPtr[0][j][i].OutputSolutionTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    for(int k=0, _dummy_parameter_=0; k<=kLastCell(); ++k)
      for(int j=0; j<=jLastCell(); ++j)
	for(int i=0; i<=iLastCell(); ++i){
	  SolnPtr[k][j][i].OutputSolutionTecplot(output_file,i,j,k);
	}
    break;
  }
}

// OutputSolutionCellTecplot
/* Outputs the solution at the centroid of every cell, including the ghost cells */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputSolutionCellTecplot(std::ofstream &output_file, const bool Title) const
{

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iLastCell() + 1 << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
    }
    
    for(int i=0, _dummy_parameter_=0; i<=iLastCell(); ++i){
      SolnPtr[0][0][i].OutputSolutionCellTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "J=" << jEnd() - jStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    }
    
    for(int j=jStart(), _dummy_parameter_=0; j<=jEnd(); ++j)
      for(int i=iStart(); i<=iEnd(); ++i){
	SolnPtr[0][j][i].OutputSolutionCellTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "J=" << jEnd() - jStart() + 1 << "\\ \n"
		  << "K=" << kEnd() - kStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
    }
    
    for(int k=kStart(); k<=kEnd(); ++k)
      for(int j=jStart(); j<=jEnd(); ++j)
	for(int i=iStart(); i<=iEnd(); ++i){
	  SolnPtr[k][j][i].OutputSolutionCellTecplot(output_file,i,j,k);
	}
    break;
  }
}

// OutputFullSolutionCellTecplot
/* Outputs the solution at the cell centroids, including the ghost cells */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputFullSolutionCellTecplot(std::ofstream &output_file, const bool Title) const
{

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << N_XYZ[0] << "\\ \n"
		  << "F = POINT \n";
    }
    
    for(int i=0, _dummy_parameter_=0; i<=iLastCell(); ++i){
	SolnPtr[0][0][i].OutputSolutionCellTecplot(output_file,i,_dummy_parameter_,_dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << N_XYZ[0] << "\\ \n"
		  << "J=" << N_XYZ[1] << "\\ \n"
		<< "F = POINT \n";
    }
    
    for(int j=0; j<=jLastCell(); ++j)
      for(int i=0; i<=iLastCell(); ++i){
	SolnPtr[0][j][i].OutputSolutionCellTecplot(output_file,i,j);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << N_XYZ[0] << "\\ \n"
		  << "J=" << N_XYZ[1] << "\\ \n"
		  << "K=" << N_XYZ[2] << "\\ \n"
		  << "F = POINT \n";
    }
    
    for(int k=0, _dummy_parameter_=0; k<=kLastCell(); ++k)
      for(int j=0; j<=jLastCell(); ++j)
	for(int i=0; i<=iLastCell(); ++i){
	  SolnPtr[k][j][i].OutputSolutionCellTecplot(output_file,i,j,k);
	}
    break;
  }
}

// OutputSolutionCellTecplotOneZone
/* Outputs the solution at the nodes without the ghost cells in a single tecplot zone */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputSolutionCellTecplotOneZone(std::ofstream &output_file, const bool Title) const
{

  int i,j,k, _dummy_parameter_(0), SubgridNy, SubgridNz;

  switch(SpaceDimension){
  case OneD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
    }

    for(i=iStart(); i<=iEnd(); ++i){
      SolnPtr[0][0][i].OutputSolutionCellTecplotOneZone(output_file, _dummy_parameter_, _dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    }

    for(j=jStart(); j<=jEnd(); ++j)
      for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
	for(i=iStart(); i<=iEnd(); ++i){
	  SolnPtr[0][j][i].OutputSolutionCellTecplotOneZone(output_file, SubgridNy, _dummy_parameter_);
	}
      }
    break;

  case ThreeD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].kSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].kSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
    }

    for(k=kStart(); k<=kEnd(); ++k)
      for(SubgridNz = 0; SubgridNz < SolnPtr[0][0][0].kSubgridPoints(); ++SubgridNz){
	for(j=jStart(); j<=jEnd(); ++j)
	  for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
	    for(i=iStart(); i<=iEnd(); ++i){
   	      SolnPtr[k][j][i].OutputSolutionCellTecplotOneZone(output_file, SubgridNy, SubgridNz);
	    }
	  }
      }
    
    break;
  }
}

// OutputPWC
/* Outputs the piecewise constant solution at the nodes of the subgrids */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
OutputPWC(std::ofstream &output_file, const bool Title) const{

  int i,j,k, _dummy_parameter_(0), SubgridNy, SubgridNz;

  switch(SpaceDimension){
  case OneD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"PWC " << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
    }

    for(i=iStart(); i<=iEnd(); ++i){
      SolnPtr[0][0][i].OutputPWC(output_file, _dummy_parameter_, _dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"PWC " << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    }

    for(j=jStart(); j<=jEnd(); ++j)
      for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
	for(i=iStart(); i<=iEnd(); ++i){
	  SolnPtr[0][j][i].OutputPWC(output_file, SubgridNy, _dummy_parameter_);
	}
      }
    break;

  case ThreeD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"PWC " << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,ThreeD);
    }

    for(k=jStart(); k<=kEnd(); ++k)
      for(SubgridNz = 0; SubgridNz < SolnPtr[0][0][0].kSubgridPoints(); ++SubgridNz){
	for(j=jStart(); j<=jEnd(); ++j)
	  for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
	    for(i=iStart(); i<=iEnd(); ++i){
	      SolnPtr[k][j][i].OutputPWC(output_file, SubgridNy, SubgridNz);
	    }
	  }
      }
    
    break;
  }
}

// OutputSolutionAtGaussPointsTecplot
/* Outputs the solution at the Gauss quadrature points in a single tecplot zone */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputSolutionAtGaussPointsTecplot(std::ofstream &output_file, const bool Title) const
{

  int i,j,k, _dummy_parameter_(0), SubgridNy, SubgridNz;

  switch(SpaceDimension){
  case OneD:
    break;

  case TwoD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"Solution " << SpaceDimension << "D" << " (at Gauss quadrature points)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*4 << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*2 << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*4 << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1)*2 << "\\ \n"
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,TwoD);
    }

    for(j=jStart(); j<=jEnd(); ++j){
      //       // Output the south Gauss points
      //       for(i=iStart(); i<=iEnd(); ++i){
      // 	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, SOUTH, _dummy_parameter_);
      //       }
      //       // Output the first west and east Gauss point
      //       for(i=iStart(); i<=iEnd(); ++i){
      // 	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, WEST, 1);
      //       }
      //       // Output the second west and east Gauss point
      //       for(i=iStart(); i<=iEnd(); ++i){
      // 	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, WEST, 2);
      //       }
      //       // Output the north Gauss point
      //       for(i=iStart(); i<=iEnd(); ++i){
      // 	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, NORTH, _dummy_parameter_);
      //       }


      // Output the row 1 of  Gauss points
      for(i=iStart(); i<=iEnd(); ++i){
      	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, 1, _dummy_parameter_);
      }

      // Output the row 2 of  Gauss points
      for(i=iStart(); i<=iEnd(); ++i){
      	SolnPtr[0][j][i].OutputSolutionAtGaussPointsTecplot(output_file, 2, _dummy_parameter_);
      }
    }
    break;

  case ThreeD:
    // --> RR: Needs to be implemented
    break;
  }
}

// OutputSmoothnessIndicatorAtNodesTecplot()
/* Outputs the smoothness indicator at the nodes of the cell */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
OutputSmoothnessIndicatorAtNodesTecplot(std::ofstream &output_file, const bool Title) const
{

  int i,j,k, _dummy_parameter_, SubgridNy, SubgridNz;

  switch(SpaceDimension){
  case OneD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"Smoothness Indicator in " << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n"
                  << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    }

    for(i=iStart(), _dummy_parameter_=0; i<=iEnd(); ++i){
      SolnPtr[0][0][i].OutputSmoothnessIndicatorAtNodesTecplot(output_file, _dummy_parameter_, _dummy_parameter_);
    }
    break;

  case TwoD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"Smoothness Indicator in " << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
                  << "\"y\" \\ \n"
                  << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    }

    for(j=jStart(), _dummy_parameter_=0; j<=jEnd(); ++j)
      for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
        for(i=iStart(); i<=iEnd(); ++i){
          SolnPtr[0][j][i].OutputSmoothnessIndicatorAtNodesTecplot(output_file, SubgridNy, _dummy_parameter_);
        }
      }
    break;

  case ThreeD:
    if (Title){
      /* Write TECPLOT title */
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
                  << "\"y\" \\ \n"
                  << "\"z\" \\ \n"
                  << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    } else {
      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1)*SolnPtr[0][0][0].iSubgridPoints() << "\\ \n"
                  << "J=" << (jEnd() - jStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "K=" << (kEnd() - kStart() + 1)*SolnPtr[0][0][0].jSubgridPoints() << "\\ \n"
                  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    }

    for(k=jStart(); k<=kEnd(); ++k)
      for(SubgridNz = 0; SubgridNz < SolnPtr[0][0][0].kSubgridPoints(); ++SubgridNz){
        for(j=jStart(); j<=jEnd(); ++j)
          for(SubgridNy = 0; SubgridNy < SolnPtr[0][0][0].jSubgridPoints(); ++SubgridNy){
            for(i=iStart(); i<=iEnd(); ++i){
              SolnPtr[k][j][i].OutputSmoothnessIndicatorAtNodesTecplot(output_file, SubgridNy, SubgridNz);
            }
          }
      }
  
    break;
  }
}

// OutputMultipleCorrelationCoefficient
/* Outputs the multiple-correlation coefficient associated with every cell center */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputMultipleCorrelationCoefficient(std::ofstream &output_file,const bool Title) const{

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "TITLE = \"Multiple-Correlation Coeffcient in " <<SpaceDimension <<"D"<< " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    }
    
    for(int i=iStart(); i<=iEnd(); ++i){
      output_file << SolnPtr[0][0][i].CellCenter() << "\t" << SolnPtr[0][0][i].CellMCC();
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "TITLE = \"Multiple-Correlation Coeffcient in " <<SpaceDimension <<"D"<< " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n"
                  << "VARIABLES = \"y\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
                  << "J=" << jEnd() - jStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \n"
                  << "VARIABLES = \"y\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
                  << "J=" << jEnd() - jStart() + 1 << "\\ \n";
      output_file << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    }
    for(int j=jStart(); j<=jEnd(); ++j)
      for(int i=iStart(); i<=iEnd(); ++i){
        output_file << SolnPtr[0][j][i].CellCenter() << "\t" << SolnPtr[0][j][i].CellMCC();
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "TITLE = \"Multiple-Correlation Coeffcient in " <<SpaceDimension <<"D"<< " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n"
                  << "VARIABLES = \"y\" \n"
                  << "VARIABLES = \"z\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
                  << "J=" << jEnd() - jStart() + 1 << "\\ \n"
                  << "K=" << kEnd() - kStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \n"
                  << "VARIABLES = \"y\" \n"
                  << "VARIABLES = \"z\" \n";
      output_file << "\"MCC\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
                  << "J=" << jEnd() - jStart() + 1 << "\\ \n"
                  << "K=" << kEnd() - kStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    }
    for(int k=kStart(); k<=kEnd(); ++k)
      for(int j=jStart(); j<=jEnd(); ++j)
        for(int i=iStart(); i<=iEnd(); ++i){
          output_file << SolnPtr[k][j][i].CellCenter() << "\t" << SolnPtr[k][j][i].CellMCC();
        }
    break;

  }
}

// OutputSmoothnessIndicatorFlag
/* Outputs the smoothness indicator flag state of every cell */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputSmoothnessIndicatorFlag(std::ofstream &output_file,const bool Title) const{
  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "TITLE = \"Smoothness Indicator Flag for " << SpaceDimension << "D" << " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n";
      output_file << "\"SI\" \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \n";
      output_file << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << iEnd() - iStart() + 1 << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE) \n";
    }
    
    for(int i=iStart(); i<=iEnd(); ++i){
      output_file << SolnPtr[0][0][i].CellCenter() << "\t" << SolnPtr[0][0][i].UnfitReconstructionFlag() << std::endl;
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "TITLE = \"Smoothness Indicator Flag for " << SpaceDimension << "D" << " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      output_file << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1) << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1) << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      output_file << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1) << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1) << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE) \n";
    }

    for(int j=jStart(); j<=jEnd(); ++j){
      for(int i=iStart(); i<=iEnd(); ++i){
	output_file << SolnPtr[0][j][i].CellCenter() << "\t" << SolnPtr[0][j][i].UnfitReconstructionFlag() << std::endl;
      }
    }
    break;
    
  case ThreeD:
    if (Title) {
      output_file << "TITLE = \"Smoothness Indicator Flag for " << SpaceDimension << "D" << " (at cell center locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
                  << "\"z\" \\ \n";
      output_file << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1) << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1) << "\\ \n"
                  << "K=" << (kEnd() - kStart() + 1) << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    }
    else {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      output_file << "\"SI\" \\ \n";
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE T = \" SolBlock \" \\ \n";
      output_file << "I=" << (iEnd() - iStart() + 1) << "\\ \n"
		  << "J=" << (jEnd() - jStart() + 1) << "\\ \n"
		  << "K=" << (kEnd() - kStart() + 1) << "\\ \n"
		  << "F = POINT \n";
      output_file << "DT = (DOUBLE DOUBLE DOUBLE DOUBLE) \n";
    }
    for(int k=kStart(); k<=kEnd(); ++k){
      for(int j=jStart(); j<=jEnd(); ++j){
        for(int i=iStart(); i<=iEnd(); ++i){
          output_file << SolnPtr[k][j][i].CellCenter() << "\t" << SolnPtr[k][j][i].UnfitReconstructionFlag() << std::endl;
        }
      }
    }
    break;
  }
}

// OutputReconstructedSolutionTecplot
/* Outputs the solution of each reconstructed function at the subgrid nodes of the domain used for its
   reconstruction */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  OutputReconstructedSolutionTecplot(std::ofstream &output_file,const bool Title) const{

  typename CompCellType::Node Node;

  switch(SpaceDimension){
  case OneD:
    if (Title) {
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    else {
      output_file << "VARIABLES = \"x\" \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }

    int cell, StencilCell, SGPoint;    
    for( cell=iStart(); cell<=iEnd(); ++cell){
      output_file << "ZONE T = \" Cell [" << cell << "] \" \\ \n"
	             /* the total number of subgrid points in the stencil used for reconstruction */
		  << "I=" << (1 + 2*SolnPtr[0][0][cell].CellRings())*SolnPtr[0][0][cell].iSubgridPoints() << "\n" 
		  << "F = POINT \n";
      VarNames.PrintDataTypeTecplot(output_file,OneD);
      for(StencilCell=cell-SolnPtr[0][0][cell].CellRings(); StencilCell<=cell+SolnPtr[0][0][cell].CellRings(); ++StencilCell){
	for(SGPoint=0; SGPoint<SolnPtr[0][0][cell].iSubgridPoints(); ++SGPoint){
	  Node = SolnPtr[0][0][StencilCell].CellSubgridSolution(SGPoint).GetNode();
	  output_file << Node << "\t"
		      << SolnPtr[0][0][cell].SolutionAtCoordinates(Node)
		      << "\n";
	}
      }
    }
    break;

  case TwoD:
    if (Title) {
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    else {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(int j=jStart(), _dummy_parameter_=0; j<=jEnd(); ++j)
      for(int i=iStart(); i<=iEnd(); ++i){
	SolnPtr[0][j][i].OutputReconstructedSolutionTecplot(output_file,i,j,_dummy_parameter_);
      }
    break;

  case ThreeD:
    if (Title) {
      output_file << "TITLE = \"Solution" << SpaceDimension << "D" << " (at subgrid nodes locations)\"\n ";
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    else {
      output_file << "VARIABLES = \"x\" \\ \n"
		  << "\"y\" \\ \n"
		  << "\"z\" \\ \n";
      VarNames.PrintHeaderTecplot(output_file);
      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);
    }
    
    for(int k=kStart(); k<=kEnd(); ++k)
      for(int j=jStart(); j<=jEnd(); ++j)
	for(int i=iStart(); i<=iEnd(); ++i){
	  SolnPtr[k][j][i].OutputReconstructedSolutionTecplot(output_file,i,j,k);
	}
    break;
  }
}

/******************************************************
void BC_Reflection(void)
Fills in the ghost cells with data solution such
that to obtain reflection at both ends of the domain
*******************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::BC_Reflection(void)
{

  switch(SpaceDimension){
  case OneD:
    for (int ghost=0; ghost<iStart(); ++ghost){
      SolnPtr[0][0][ghost].CellSolution() = -SolnPtr[0][0][2*iStart()-ghost-1].CellSolution();
      SolnPtr[0][0][iEnd()+ghost+1].CellSolution() = -1.0*SolnPtr[0][0][iEnd()-ghost].CellSolution();

      /* Update the first derivative */
      SolnPtr[0][0][ghost].CellDeriv(0) = SolnPtr[0][0][ghost].CellSolution();
      SolnPtr[0][0][iEnd()+ghost+1].CellDeriv(0) = SolnPtr[0][0][iEnd()+ghost+1].CellSolution();
    }
    UpdateSubgridSolution();
    break;
  case TwoD:
    break;

  case ThreeD:
    break;
  }

}

/***********************************************
void ReconstructSolution(Reconstruct1D_Input_Parameters)
Solves the reconstruction over the domain
***********************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructSolution(const Reconstruct1D_Input_Parameters & IP)
{
  int i;
  //******* For piecewise constant ********
  if(SolnPtr[0][0][0].CellRings() == 0){
    for (i=iStart(); i<=iEnd(); ++i){
      SolnPtr[0][0][i].CellDeriv()(0) = SolnPtr[0][0][i].CellSolution();
      SolnPtr[0][0][i].UpdateSubgridSolution();
    }
    return;
  }

  vector<int> i_index(1 + 2*SolnPtr[0][0][0].CellRings()); // CellsInOneDirection = 1 + 2*Soln[0].CellRings();
  int UnfitCells(0); 		/* the number of cells that didn't get smooth reconstruction */
  int TotalModifiedCells(0); 	/* the total number of cells affected by the order reduction */

  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/
 
  switch(IP.ReconstructionMethod()){
    
    /* DataDependent ENO-like reconstruction */
  case DD_ENO:

    break;
    
    /* Essentially NonOscillatory (ENO) */
  case ENO:
    for (i=iStart(); i<=iEnd(); ++i){
      /* Solve reconstruction for the current cell */
      ENO_Reconstruction(*this,i);
    }
    break;

    /* Essentially NonOscillatory (ENO) with Least-Squares */
  case ENO_LS:
    std::cout << "\n This reconstruction type is not used anymore.\n"
	      << " Still interested in this reconstruction type?\n"
	      << " Check out versions of this file up to 2.20.\n";
    exit(1);
    break;

    /* Weighted Essentially NonOscillatory (WENO) */
  case WENO:
    /* Haselbacher's reconstruction */
    std::cout << "\n This reconstruction type is not used anymore.\n"
	      << " Still interested in this reconstruction type?\n"
	      << " Check out versions of this file up to 2.20.\n";
    exit(1);
    break;

    /* Spectral Differences */
  case SpectralDiff:
    
    break;

  case CENO:
    /************* STEP 1: Reconstruction the solution in every cell ***************/
    for (i=iStart()-SolnPtr[0][0][0].CellRings(); i<=iEnd()+SolnPtr[0][0][0].CellRings(); ++i){
      /* Make Stencil */
      MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,i_index);
      /* Solve reconstruction for the current cell */
      kExact_Reconstruction(*this,i_index,i,NumberOfTaylorDerivatives(),OFF);
    }

    /************* STEP 2: Compute the Multiple-Correlation Coefficient ************/
    for (i=iStart(); i<=iEnd(); ++i){
      /* Make Stencil */
      MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,i_index);
      
      /* Estimate the correlation coefficient for the current cell */
      ComputeSmoothnessIndicator(*this,i_index,i);
    }

    /************* STEP 3: Flag the cells which didn't get a good reconstruction ***********/
    for (i=iStart(); i<=iEnd(); ++i){
      if(SolnPtr[0][0][i].CellMCC() < CENO_Tolerances::Fit_Tolerance ){

	/* Flag the cell with non-smooth reconstruction */
	SolnPtr[0][0][i].UnfitReconstructionFlag() = ON;

#ifdef CENO_Padding
	/* flag one cell to the right and one cell to the left */
	SolnPtr[0][0][i-1].UnfitReconstructionFlag() = ON;
	SolnPtr[0][0][i+1].UnfitReconstructionFlag() = ON;
#endif
	/* Update the number of cell with not enough resolution */
	++UnfitCells;
      }
    }

    /************* STEP 4: Reduce to limited piecewise linear the reconstruction of the flagged cells ***********/
    for (i=iStart(); i<=iEnd(); ++i){
      if(SolnPtr[0][0][i].UnfitReconstructionFlag() == ON ){
	++TotalModifiedCells;
	ReconstructToLimitedPiecewiseLinear(i,IP.Limiter());
      }
    }

    /* Print UnfitCells and TotalModifiedCells */
    std::cout << std::endl;
    Print_(UnfitCells);
    Print_(TotalModifiedCells);
    break;

  default:
    std::cout << "\n Unknown high-order reconstruction.\n "
	      << "Please check in the input file the \"Reconstruction_Type\" field !\n";
    exit(1);
  }
}


/***********************************************
void ReconstructToLimitedPiecewiseLinear(iCell)
Solves the reconstruction over the domain
***********************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructToLimitedPiecewiseLinear(int iCell,
													  const int Limiter){

  vector<int> i_index(3);
  int i, n, n2, n_pts, index[2];
  double u0Min, u0Max, uQuad[2];

  n_pts = 2;
  index[0] = iCell-1;
  index[1] = iCell+1; 


  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/

  /* Make Stencil */
  MakeReconstructionStencil(1,iCell,i_index);

  /* Solve the PWL reconstruction for the current cell */
  kExact_Reconstruction(*this,i_index,iCell,2,ON);

  /* Copy the derivatives in the TD container */
  SolnPtr[0][0][iCell].CellDeriv(0,true,true,true).D() = SolnPtr[0][0][iCell].CellDerivFirstOrder(0,true,true,true).D();
  SolnPtr[0][0][iCell].CellDeriv(1,true,true,true).D() = SolnPtr[0][0][iCell].CellDerivFirstOrder(1,true,true,true).D();
  for (int TD = 2; TD<SolnPtr[0][0][iCell].NumberOfTaylorDerivatives(); ++TD){
    SolnPtr[0][0][iCell].CellDeriv(TD,true,true,true).D() = SolutionType(0.0);
  }

  /* Compute and assign the limiter */
  u0Min = SolnPtr[0][0][iCell].CellSolution();
  u0Max = u0Min;
  for ( n2 = 0 ; n2 <= 1 ; ++n2 ) {
    u0Min = min(u0Min, SolnPtr[0][0][ index[n2] ].CellSolution());
    u0Max = max(u0Max, SolnPtr[0][0][ index[n2] ].CellSolution());
  } /* endfor */

  uQuad[0] = SolnPtr[0][0][iCell].SolutionAt(-HALF*SolnPtr[0][0][iCell].CellDelta());
  uQuad[1] = SolnPtr[0][0][iCell].SolutionAt(HALF*SolnPtr[0][0][iCell].CellDelta());

  switch(Limiter) {
  case LIMITER_ZERO :
    SolnPtr[0][0][iCell].CellLimiter() = ZERO;
    break;
  case LIMITER_ONE :
    SolnPtr[0][0][iCell].CellLimiter() = ONE;
    break;
  case LIMITER_BARTH_JESPERSEN :
    SolnPtr[0][0][iCell].CellLimiter() = Limiter_BarthJespersen(uQuad, SolnPtr[0][0][iCell].CellSolution(), u0Min, u0Max, 2);
    break;
  case LIMITER_VENKATAKRISHNAN :
    SolnPtr[0][0][iCell].CellLimiter() = Limiter_Venkatakrishnan(uQuad, SolnPtr[0][0][iCell].CellSolution(), u0Min, u0Max, 2);
    break;
  case LIMITER_VANLEER :
    SolnPtr[0][0][iCell].CellLimiter() = Limiter_VanLeer(uQuad, SolnPtr[0][0][iCell].CellSolution(), u0Min, u0Max, 2);
    break;
  case LIMITER_VANALBADA :
    SolnPtr[0][0][iCell].CellLimiter() = Limiter_VanAlbada(uQuad, SolnPtr[0][0][iCell].CellSolution(), u0Min, u0Max, 2);
    break;
  default:
    SolnPtr[0][0][iCell].CellLimiter() = Limiter_BarthJespersen(uQuad, SolnPtr[0][0][iCell].CellSolution(), u0Min, u0Max, 2);
    break;
  } /* endswitch */

}

/***********************************************
void Reconstruct2DSolution(Reconstruct2D_Input_Parameters)
Solves the reconstruction over the domain
***********************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructSolution(const Reconstruct2D_Input_Parameters  & IP)
{

  require(SpaceDimension==TwoD, "ERROR:Reconstruct2DSolution(). This function can be called only for a 2D computational domain. Check your space dimension and call for the proper function!");

  int ProgressFrequency = ((iEnd()-iStart()+1)*(jEnd()-jStart()+1))/10;
  int CellsWithSmallerOrder = 0;

  if(SolnPtr[0][0][0].CellRings() == 0){
    // Piecewise constant
    for (int i=iStart(); i<=iEnd(); ++i)
      for (int j=jStart(); j<=jEnd(); ++j){
	SolnPtr[0][j][i].CellDeriv()(0,0) = SolnPtr[0][j][i].CellSolution();
	SolnPtr[0][j][i].UpdateSubgridSolution();
      }
    return;
  }

  int *i_index(NULL), *j_index(NULL);
  int UnfitCells(0); 		/* the number of cells that didn't get the smooth resolution */
  int TotalModifiedCells(0); 	/* the total number of cells affected by the order reduction */

  int i,j, FinishedCell;

  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/

  /* Determine the number of cells in the stencil based on the number of rings */
  int NumOfCellsInOneDirection, StencilSize;
  
  switch(IP.ReconstructionMethod()){

  case CENO:
    NumOfCellsInOneDirection = 2*SolnPtr[0][0][0].CellRings() + 1;
    StencilSize = NumOfCellsInOneDirection*NumOfCellsInOneDirection;
    i_index = new int [StencilSize];
    j_index = new int [StencilSize];

    // Solve reconstruction for each cell in the domain plus 2 additional layers of cells for each boundary.
    // These boundary cells are used for checking the goodness of fit of the first domain cell.
    for (i=iStart()-2, FinishedCell=0; i<=iEnd()+2; ++i){
      for (j=jStart()-2; j<=jEnd()+2; ++j){
	/* Make Stencil */
	MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,i_index,j_index);
	/* Solve reconstruction for the current cell */
	kExact_Reconstruction(*this,i_index,j_index,StencilSize);

	++FinishedCell;
	Print_Progress(FinishedCell,ProgressFrequency);
      }//endfor(j)
    }//endfor(i)


    // Compute the Smoothness Indicator for each cell reconstruction
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j){
	/* Make Stencil */
	MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,i_index,j_index);
	ComputeSmoothnessIndicator(*this,i_index,j_index,StencilSize,i,j);
    }

    // Flag the cells which didn't get a good reconstruction
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j){
	if ( SolnPtr[0][j][i].CellMCC() < CENO_Tolerances::Fit_Tolerance ){

	  /* Flag the cell with non-smooth reconstruction */
	  SolnPtr[0][j][i].UnfitReconstructionFlag() = ON;

#ifdef CENO_Padding
	  /* Flag all the cells surrounding this cell */
	  SolnPtr[0][j-1][i-1].UnfitReconstructionFlag() = ON;
	  SolnPtr[0][j-1][i].UnfitReconstructionFlag() = ON;
	  SolnPtr[0][j-1][i+1].UnfitReconstructionFlag() = ON;

	  SolnPtr[0][j][i-1].UnfitReconstructionFlag() = ON;
	  SolnPtr[0][j][i+1].UnfitReconstructionFlag() = ON;

	  SolnPtr[0][j+1][i-1].UnfitReconstructionFlag() = ON;
	  SolnPtr[0][j+1][i].UnfitReconstructionFlag() = ON;
	  SolnPtr[0][j+1][i+1].UnfitReconstructionFlag() = ON;
#endif

	  /* Update the number of cell with not enough resolution */
	  ++UnfitCells;
	}
      }//endfor

    /* reduce to limited piecewise linear */
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j){
	if(SolnPtr[0][j][i].UnfitReconstructionFlag() == ON ){
	  ++TotalModifiedCells;
	  ReconstructToLimitedPiecewiseLinear(i,j,IP.Limiter());
	}
      }

    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j){
	// Update subgrid solution
	SolnPtr[0][j][i].UpdateSubgridSolution();
      }

    /* Print UnfitCells and TotalModifiedCells */
    std::cout << std::endl;
    Print_(UnfitCells);
    Print_(TotalModifiedCells);
    break;
   
    /* DataDependent ENO-like reconstruction */
  case DD_ENO:

    break;

    /* Essentially NonOscillatory (ENO) */
  case ENO:

    break;

    /* Weighted Essentially NonOscillatory (WENO) */
  case WENO:

    break;

    /* Spectral Differences */
  case SpectralDiff:

    break;

  default:
    std::cout << "\n Unknown high-order reconstruction.\n "
	      << "Please check in the input file the \"Reconstruction_Type\" field !\n";
    exit(1);
  }

  // Free memory
  delete [] i_index; i_index = NULL;
  delete [] j_index; j_index = NULL;
}

/***********************************************************************************************
void Reconstruct2DZonalSolution(Reconstruct2D_Input_Parameters)
Solves the reconstruction over a restricted domain (used for testing external reconstructions)
***********************************************************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructZonalSolution(const Reconstruct2D_Input_Parameters  & IP)
{

  int ProgressFrequency = ((iEnd()-iStart()+1)*(jEnd()-jStart()+1))/10;

  if(SolnPtr[0][0][0].CellRings() == 0){
    // Piecewise constant
    for (int i=iStart(); i<=iEnd(); ++i)
      for (int j=jStart(); j<=jEnd(); ++j){
	SolnPtr[0][j][i].CellDeriv()(0,0) = SolnPtr[0][j][i].CellSolution();
	SolnPtr[0][j][i].UpdateSubgridSolution();
      }
    return;
  }

  int *i_index(NULL), *j_index(NULL);
  int UnfitCells(0); 		/* the number of cells that didn't get the smooth resolution */
  int TotalModifiedCells(0); 	/* the total number of cells affected by the order reduction */

  int i,j, FinishedCell;

  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/

  /* Determine the number of cells in the stencil based on the number of rings */
  int NumOfCellsInOneDirection, StencilSize;
  
  NumOfCellsInOneDirection = 2*SolnPtr[0][0][0].CellRings() + 1;
  StencilSize = NumOfCellsInOneDirection*NumOfCellsInOneDirection;
  i_index = new int [StencilSize];
  j_index = new int [StencilSize];

  // Solve reconstruction for each cell in the domain
  for (i=iStart(), FinishedCell=0; i<=iEnd(); ++i){
    for (j=jStart(); j<=jEnd(); ++j){
      /* Make Stencil */
      MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,i_index,j_index);
      /* Solve reconstruction for the current cell */
      kExact_Reconstruction(*this,i_index,j_index,StencilSize);

      Print_3(i,j,SolnPtr[0][j][i].CellDeriv());

      ++FinishedCell;
      Print_Progress(FinishedCell,ProgressFrequency);
    }//endfor(j)
  }//endfor(i)

  // Compute the Smoothness Indicator for each cell reconstruction
  for (i=iStart()+2; i<=iEnd()-2; ++i)
    for (j=jStart()+2; j<=jEnd()-2; ++j){
      /* Make Stencil */
      MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,i_index,j_index);
      ComputeSmoothnessIndicator(*this,i_index,j_index,StencilSize,i,j);
      std::cout << std::endl;
      Print_(SolnPtr[0][j][i].CellMCC());
    }

  Print_(CENO_Tolerances::Fit_Tolerance)

  // Flag the cells which didn't get a good reconstruction
  for (i=iStart()+2; i<=iEnd()-2; ++i)
    for (j=jStart()+2; j<=jEnd()-2; ++j){
      if ( SolnPtr[0][j][i].CellMCC() < CENO_Tolerances::Fit_Tolerance ){
	SolnPtr[0][j][i].UnfitReconstructionFlag() = ON;
	/* Update the number of cells with not enough resolution */
	++UnfitCells;
      }
    }//endfor

  /* reduce to limited piecewise linear */
  for (i=iStart(); i<=iEnd(); ++i)
    for (j=jStart(); j<=jEnd(); ++j){
      if(SolnPtr[0][j][i].UnfitReconstructionFlag() == ON ){
	++TotalModifiedCells;
	ReconstructToLimitedPiecewiseLinear(i,j,IP.Limiter());
      }
    }

  for (i=iStart(); i<=iEnd(); ++i)
    for (j=jStart(); j<=jEnd(); ++j){
      // Update subgrid solution
      SolnPtr[0][j][i].UpdateSubgridSolution();
    }

  /* Print UnfitCells and TotalModifiedCells */
  std::cout << std::endl;
  Print_(UnfitCells);
  Print_(TotalModifiedCells);
   
  // Free memory
  delete [] i_index; i_index = NULL;
  delete [] j_index; j_index = NULL;
}

/*****************************************************
void ReconstructToLimitedPiecewiseLinear(iCell,jCell)
Solves the reconstruction over the domain
******************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructToLimitedPiecewiseLinear(int iCell, int jCell,
													  const int Limiter)
{

#ifdef LimiterBasedOnlyOnGaussPoints
  int NQuad(8);
  double uQuad[8];
#else
  int NQuad(12);
  double uQuad[12];
#endif

  int i, n, n2, n_pts, indexI[8], indexJ[8];
  double u0Min, u0Max;

  n_pts = 8;
  /* First ring */
  indexI[0] = iCell-1; indexJ[0] = jCell-1;
  indexI[1] = iCell  ; indexJ[1] = jCell-1;
  indexI[2] = iCell+1; indexJ[2] = jCell-1;
  indexI[3] = iCell-1; indexJ[3] = jCell;
  indexI[4] = iCell+1; indexJ[4] = jCell;
  indexI[5] = iCell-1; indexJ[5] = jCell+1;
  indexI[6] = iCell  ; indexJ[6] = jCell+1;
  indexI[7] = iCell+1; indexJ[7] = jCell+1;
 
  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/
  int *i_index, *j_index, StencilSize(9);
  i_index = new int [StencilSize];
  j_index = new int [StencilSize];

  /* Make Stencil */
  MakeReconstructionStencil(1,iCell,jCell,i_index,j_index);

  /* Zero all the high-order derivatives */
  for (i=1; i<SolnPtr[0][jCell][iCell].NumberOfTaylorDerivatives(); ++i){
    SolnPtr[0][jCell][iCell].CellDeriv(i,true,true,true).D() = SolutionType(0.0);
  }

  /* Solve reconstruction for the current cell */
  FirstOrder_kExact_Reconstruction(*this, i_index, j_index, StencilSize);

  /* Compute and assign the limiter */
  u0Min = SolnPtr[0][jCell][iCell].CellSolution();
  u0Max = u0Min;
  for ( n2 = 0 ; n2 < n_pts ; ++n2 ) {
    u0Min = min(u0Min, SolnPtr[0][ indexJ[n2] ][ indexI[n2] ].CellSolution());
    u0Max = max(u0Max, SolnPtr[0][ indexJ[n2] ][ indexI[n2] ].CellSolution());
  } /* endfor */

  /* Gauss Quadrature Point values */
  Vector2D GP1, GP2;

  /* N face */
  SolnPtr[0][jCell][iCell].CellGeometry().GaussQuadPointsFaceN(GP1,GP2);
  uQuad[0] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);
  uQuad[1] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP2.x, GP2.y);

  /* S face */
  SolnPtr[0][jCell][iCell].CellGeometry().GaussQuadPointsFaceS(GP1,GP2);
  uQuad[2] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);
  uQuad[3] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP2.x, GP2.y);

  /* E face */
  SolnPtr[0][jCell][iCell].CellGeometry().GaussQuadPointsFaceE(GP1,GP2);
  uQuad[4] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);
  uQuad[5] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP2.x, GP2.y);

  /* W face */
  SolnPtr[0][jCell][iCell].CellGeometry().GaussQuadPointsFaceW(GP1,GP2);
  uQuad[6] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);
  uQuad[7] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP2.x, GP2.y);

#ifndef LimiterBasedOnlyOnGaussPoints
  /* Use also the nodes of the cell */
  GP1 = SolnPtr[0][jCell][iCell].CellGeometry().NodeNW();
  uQuad[8] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);

  GP1 = SolnPtr[0][jCell][iCell].CellGeometry().NodeNE();
  uQuad[9] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);

  GP1 = SolnPtr[0][jCell][iCell].CellGeometry().NodeSW();
  uQuad[10] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);

  GP1 = SolnPtr[0][jCell][iCell].CellGeometry().NodeSE();
  uQuad[11] = SolnPtr[0][jCell][iCell].SolutionAtCoordinates(GP1.x, GP1.y);
#endif

  switch(Limiter) {
  case LIMITER_ZERO :
    SolnPtr[0][jCell][iCell].CellLimiter() = ZERO;
    break;
  case LIMITER_ONE :
    SolnPtr[0][jCell][iCell].CellLimiter() = ONE;
    break;
  case LIMITER_BARTH_JESPERSEN :
    SolnPtr[0][jCell][iCell].CellLimiter()=Limiter_BarthJespersen(uQuad,SolnPtr[0][jCell][iCell].CellSolution(),u0Min,u0Max,NQuad);
    break;
  case LIMITER_VENKATAKRISHNAN :
    SolnPtr[0][jCell][iCell].CellLimiter() = Limiter_Venkatakrishnan(uQuad,SolnPtr[0][jCell][iCell].CellSolution(),u0Min,u0Max,NQuad);
    break;
  case LIMITER_VANLEER :
    SolnPtr[0][jCell][iCell].CellLimiter() = Limiter_VanLeer(uQuad,SolnPtr[0][jCell][iCell].CellSolution(),u0Min,u0Max,NQuad);
    break;
  case LIMITER_VANALBADA :
    SolnPtr[0][jCell][iCell].CellLimiter() = Limiter_VanAlbada(uQuad,SolnPtr[0][jCell][iCell].CellSolution(),u0Min,u0Max,NQuad);
    break;
  default:
    SolnPtr[0][jCell][iCell].CellLimiter() = Limiter_BarthJespersen(uQuad,SolnPtr[0][jCell][iCell].CellSolution(),u0Min,u0Max,NQuad);
    break;
  } /* endswitch */

  // Free memory
  delete [] i_index; i_index = NULL;
  delete [] j_index; j_index = NULL;
}

/***********************************************
void ReconstructSolution(Reconstruct3D_Input_Parameters)
Solves the reconstruction over the domain
***********************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructSolution(const Reconstruct3D_Input_Parameters & IP)
{

  require(SpaceDimension==ThreeD, "ERROR:Reconstruct3DSolution(). This function can be called only for a 3D computational domain. Check your space dimension and call for the proper function!");

  int ProgressFrequency = (iEnd()-iStart()+1)*(jEnd()-jStart()+1)*(kEnd()-kStart()+1)/10;
  int CellsWithSmallerOrder = 0;

  if(SolnPtr[0][0][0].CellRings() == 0){
    // Piecewise constant
    for (int i=iStart(); i<=iEnd(); ++i)
      for (int j=jStart(); j<=jEnd(); ++j)
	for (int k=kStart(); k<=kEnd(); ++k){
	  SolnPtr[k][j][i].CellDeriv()(0,0,0) = SolnPtr[k][j][i].CellSolution();
	  SolnPtr[k][j][i].UpdateSubgridSolution();
	}
    return;
  }

  //vector<int> i_index, j_index, k_index;
  int *i_index(NULL), *j_index(NULL), *k_index(NULL);
  int UnfitCells(0);            /* the number of cells that didn't get the smooth resolution */
  int TotalModifiedCells(0);    /* the total number of cells affected by the order reduction */

  int i,j,k, FinishedCell;

  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/

  /* Determine the number of cells in the stencil based on the number of rings */
  int NumOfCellsInOneDirection, StencilSize;
  
  switch(IP.ReconstructionMethod()){
    
  case CENO:
    NumOfCellsInOneDirection = 2*SolnPtr[0][0][0].CellRings() + 1;
    StencilSize = NumOfCellsInOneDirection*NumOfCellsInOneDirection*NumOfCellsInOneDirection;
    i_index = new int [StencilSize];
    j_index = new int [StencilSize];
    k_index = new int [StencilSize];
    
    // Solve reconstruction for each cell in the domain plus 2 additional layers of cells for each boundary.
    // These boundary cells are used for checking the goodness of fit of the first domain cell.
    for (i=iStart()-2, FinishedCell=0; i<=iEnd()+2; ++i){
      for (j=jStart()-2; j<=jEnd()+2; ++j){
	for (k=kStart()-2; k<=kEnd()+2; ++k){
	  /* Make Stencil */
	  MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,k,i_index,j_index,k_index);
	  /* Solve reconstruction for the current cell */
	  kExact_Reconstruction(*this,i_index,j_index,k_index,StencilSize);
	  ++FinishedCell;
	  Print_Progress(FinishedCell,ProgressFrequency);
	}//endfor(k)
      }//endfor(j)
    }//endfor(i)
    
    
    // Compute the Smoothness Indicator for each cell reconstruction

 
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j)
	for (k=kStart(); k<=kEnd(); ++k){
	  /* Make Stencil */
	  MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,j,k,i_index,j_index,k_index);
	  ComputeSmoothnessIndicator(*this,i_index,j_index,k_index,StencilSize,i,j,k); 
	}

    // Flag the cells which didn't get a good reconstruction
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j)
	for (k=kStart(); k<=kEnd(); ++k){
	  if ( SolnPtr[k][j][i].CellMCC() < CENO_Tolerances::Fit_Tolerance ){
            /* Flag the cell with non-smooth reconstruction */
            SolnPtr[k][j][i].UnfitReconstructionFlag() = ON;
	  
            //#ifdef CENO_Padding
            //	  /* Flag all the cells surrounding this cell */
            //	  SolnPtr[0][j-1][i-1].UnfitReconstructionFlag() = ON;
            //	  SolnPtr[0][j-1][i].UnfitReconstructionFlag() = ON;
            //	  SolnPtr[0][j-1][i+1].UnfitReconstructionFlag() = ON;
            //
            //	  SolnPtr[0][j][i-1].UnfitReconstructionFlag() = ON;
            //	  SolnPtr[0][j][i+1].UnfitReconstructionFlag() = ON;
            //
            //	  SolnPtr[0][j+1][i-1].UnfitReconstructionFlag() = ON;
            //	  SolnPtr[0][j+1][i].UnfitReconstructionFlag() = ON;
            //	  SolnPtr[0][j+1][i+1].UnfitReconstructionFlag() = ON;
            //#endif

            /* Update the number of cell with not enough resolution */
            ++UnfitCells;
          }//end if
        }//endfor
    //
//    /* reduce to limited piecewise linear */
//    for (i=iStart(); i<=iEnd(); ++i)
//      for (j=jStart(); j<=jEnd(); ++j)
//	for (k=kStart(); k<=kEnd(); ++k){
//	  if(SolnPtr[k][j][i].UnfitReconstructionFlag() == ON ){
//	    ++TotalModifiedCells;
//	    ReconstructToLimitedPiecewiseLinear(i,j,k,IP.Limiter());  // --> RR: Must be Created
//	  }
//	}
//
    for (i=iStart(); i<=iEnd(); ++i)
      for (j=jStart(); j<=jEnd(); ++j)
	for (k=kStart(); k<=kEnd(); ++k){
	// Update subgrid solution
	SolnPtr[k][j][i].UpdateSubgridSolution();
        }

    /* Print UnfitCells and TotalModifiedCells */
    std::cout << std::endl;
    Print_(UnfitCells);
    Print_(TotalModifiedCells);

    // Free memory
    delete [] i_index; i_index = NULL;
    delete [] j_index; j_index = NULL;
    delete [] k_index; k_index = NULL;
    break;

    /* DataDependent ENO-like reconstruction */
  case DD_ENO:
    
    break;

    /* Essentially NonOscillatory (ENO) */
  case ENO:

    break;

    /* Weighted Essentially NonOscillatory (WENO) */
  case WENO:

    break;

    /* Spectral Differences */
  case SpectralDiff:

    break;
  }

}

/***********************************************************************************
void ReconstructSolutionAndBoundary(Reconstruct1D_Input_Parameters)
Solves the reconstruction over the domain and in the first ghost cell at every end
************************************************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
  void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::ReconstructSolutionAndBoundary(const Reconstruct1D_Input_Parameters & IP)
{

  require(SpaceDimension==OneD, "ERROR:Reconstruct1DSolutionAndBoundary(). This function can be called only for a 1D computational domain. Check your space dimension and call for the proper function!");

  //  int ProgressFrequency = (iEnd()-iStart()+1)/10;
  //  int CellsWithSmallerOrder = 0;

  if(SolnPtr[0][0][0].CellRings() == 0){
    // Piecewise constant
    for (int i=iStart()-1; i<=iEnd()+1; ++i){
      SolnPtr[0][0][i].CellDeriv()(0) = SolnPtr[0][0][i].CellSolution();
      SolnPtr[0][0][i].UpdateSubgridSolution();
    }
    return;
  }

  /* Check if there are enough ghost cells */
  require(iStart()-SolnPtr[0][0][0].CellRings()-1 >= 0, "ERROR:Reconstruct1DSolutionAndBoundary(). Not enough ghost cells!!");

  vector<int> i_index;

  /*****************************************************************************
   *                    COMPUTE THE RECONSTRUCTION                             *
   *****************************************************************************/

  /* Determine the number of cells in the stencil based on the number of rings */
  int CellsInOneDirection;
  
  switch(IP.ReconstructionMethod()){
    
    /* DataDependent ENO-like reconstruction */
  case DD_ENO:

    break;
    
    /* Essentially NonOscillatory (ENO) */
  case ENO:
    
    for (int i=iStart()-1; i<=iEnd()+1; ++i){
      /* Solve reconstruction for the current cell */
      ENO_Reconstruction(*this,i);
    }
    break;

    /* Essentially NonOscillatory (ENO) with Least-Squares */
  case ENO_LS:

    break;

    /* Weighted Essentially NonOscillatory (WENO) */
  case WENO:
    
    break;

    /* Spectral Differences */
  case SpectralDiff:
    
    break;

  default:
    std::cout << "\n Unknown high-order reconstruction.\n "
	      << "Please check in the input file the \"Reconstruction_Type\" field !\n";
    exit(1);
  }

  // Update subgrid solution
  UpdateSubgridSolution();
}

/****************************************************************
SetInitialData()
****************************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
  template<class InputParameters> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetInitialData(const InputParameters & IP)
{

  int ProgressFrequency;

  switch(SpaceDimension){
  case OneD:
    ProgressFrequency = IP.iCell()/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=0, FinishedCell=0; i<=iLastCell(); ++i){
      SolnPtr[0][0][i].SetInitialSolution(IP.TestFunction());
      /* Update subgrid solution */
      SolnPtr[0][0][i].UpdateSubgridSolution();
      /* Print progress */
      ++FinishedCell;
      Print_Progress(FinishedCell,ProgressFrequency);
    }
    std::cout << "\n";
    /* Determine the MaxDeltaSolutionOverDomain */
    SetMaxDeltaSolutionOverDomain();
    break;
  case TwoD:
    ProgressFrequency  = (IP.iCell()*IP.jCell())/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=0, FinishedCell=0; i<=iLastCell(); ++i)
      for (int j=0; j<=jLastCell(); ++j){
	SolnPtr[0][j][i].SetInitialSolution(IP.TestFunction());
	/* Update subgrid solution */
	SolnPtr[0][j][i].UpdateSubgridSolution();
	/* Print progress */
	++FinishedCell;
	Print_Progress(FinishedCell,ProgressFrequency);
      }
    std::cout << "\n";
    /* Determine the MaxDeltaSolutionOverDomain */
    SetMaxDeltaSolutionOverDomain();
    break;
  case ThreeD:
    ProgressFrequency  = (IP.iCell()*IP.jCell()*IP.kCell())/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=0, FinishedCell=0; i<=iLastCell(); ++i)
      for (int j=0; j<=jLastCell(); ++j)
	for (int k=0; k<=kLastCell(); ++k){
	  SolnPtr[k][j][i].SetInitialSolution(IP.TestFunction());
	  /* Update subgrid solution */
	  SolnPtr[k][j][i].UpdateSubgridSolution();
	  /* Print progress */
	  ++FinishedCell;
	  Print_Progress(FinishedCell,ProgressFrequency); 
	}
    std::cout << "\n";
    /* Determine the MaxDeltaSolutionOverDomain */
    SetMaxDeltaSolutionOverDomain();
    break;
  }

}

// UpdateSubgridSolution(void)
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  UpdateSubgridSolution(void)
{
  switch(SpaceDimension){
  case OneD:
    for (int i=0; i<=iLastCell(); ++i){
      SolnPtr[0][0][i].UpdateSubgridSolution();
    }
    break;
  case TwoD:
    for (int i=0; i<=iLastCell(); ++i)
      for (int j=0; j<=jLastCell(); ++j){
	SolnPtr[0][j][i].UpdateSubgridSolution();
      }
    break;
  case ThreeD:
    for (int i=0; i<=iLastCell(); ++i)
      for (int j=0; j<=jLastCell(); ++j)
	for (int k=0; k<=kLastCell(); ++k){
	  SolnPtr[k][j][i].UpdateSubgridSolution();
	}
    break;
  }
}

// SetMaxDeltaSolutionOverDomain(void)
/* Computes the difference between MaxSolution and MinSolution over the whole domain 
   The MaxDeltaSolutionOverDomain parameter is used in the computation of DD weights during
   the reconstruction.*/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  SetMaxDeltaSolutionOverDomain(void)
{

  ColumnVector MinSolution(NumberOfParameters), MaxSolution(NumberOfParameters);

  switch(SpaceDimension){
  case OneD:

    for(int parameter=1; parameter <= NumberOfParameters; ++parameter){
      /* Initialize MaxSolution and MinSolution for every parameter */
      MaxSolution(parameter-1) = SolnPtr[0][0][0].CellSolution(parameter);
      MinSolution(parameter-1) = SolnPtr[0][0][0].CellSolution(parameter);
      
      /* Determine the MinSolution and MaxSolution */
      for (int i=0; i<=iLastCell(); ++i){ // for each cell of the computational domain
	MinSolution(parameter-1) = min(MinSolution(parameter-1),SolnPtr[0][0][i].CellSolution(parameter));
	MaxSolution(parameter-1) = max(MaxSolution(parameter-1),SolnPtr[0][0][i].CellSolution(parameter));
      }
      /* Compute MaxDeltaSolutionOverDomain */
      MaxDeltaSolutionOverDomain[parameter-1] = MaxSolution(parameter-1) - MinSolution(parameter-1);
    }

    break;
  case TwoD:

    for(int parameter=1; parameter <= NumberOfParameters; ++parameter){
      /* Initialize MaxSolution and MinSolution for every parameter */
      MaxSolution(parameter-1) = SolnPtr[0][0][0].CellSolution(parameter);
      MinSolution(parameter-1) = SolnPtr[0][0][0].CellSolution(parameter);
    
      /* Determine the MinSolution and MaxSolution for every parameter */
      for (int i=0; i<=iLastCell(); ++i)
	for (int j=0; j<=jLastCell(); ++j){ // for each cell of the computational domain
	  MinSolution(parameter-1) = min(MinSolution(parameter-1),SolnPtr[0][j][i].CellSolution(parameter));
	  MaxSolution(parameter-1) = max(MaxSolution(parameter-1),SolnPtr[0][j][i].CellSolution(parameter));
	}

      /* Compute MaxDeltaSolutionOverDomain for every parameter */
      MaxDeltaSolutionOverDomain[parameter-1] = MaxSolution(parameter-1) - MinSolution(parameter-1);
    }

    break;
  case ThreeD:
    for(int parameter=1; parameter <= NumberOfParameters; ++parameter){
      /* Initialize MaxSolution and MinSolution */
      MaxSolution(parameter-1) = SolnPtr[kStart()][jStart()][iStart()].CellSolution(parameter);
      MinSolution(parameter-1) = SolnPtr[kStart()][jStart()][iStart()].CellSolution(parameter);

      /* Determine the MinSolution and MaxSolution */
      for (int i=0; i<=iLastCell(); ++i) 
	for (int j=0; j<=jLastCell(); ++j)
	  for (int k=0; k<=kLastCell(); ++k){ // for each cell of the computational domain
	    MinSolution(parameter-1) = min(MinSolution(parameter-1),SolnPtr[k][j][i].CellSolution(parameter));
	    MaxSolution(parameter-1) = max(MaxSolution(parameter-1),SolnPtr[k][j][i].CellSolution(parameter));
	  }
    
      /* Compute MaxDeltaSolutionOverDomain */
      MaxDeltaSolutionOverDomain[parameter-1] = MaxSolution(parameter-1) - MinSolution(parameter-1);
    }

    break;
  }
}

/********************************************************
* template<class InputParameters>                       *
* AssessReconstructionAccuracy(const InputParameters)   *
*                                                       *
* Computes the reconstruction error in each cell and    *
* then computes the L1_norm, L2_norm and L_max.         *
*********************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
  template<class InputParameters> inline
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  AssessReconstructionAccuracy(const InputParameters & IP)
{
  int ProgressFrequency;
  double Domain(0.0);

  switch(SpaceDimension){
  case OneD:
    ProgressFrequency = IP.iCell()/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=iStart(), FinishedCell=0; i<=iEnd(); ++i){
      /* compute error in each cell */
      SolnPtr[0][0][i].ComputeReconstructionError(IP.TestFunction());
      /* update norms */
      L1Norm += SolnPtr[0][0][i].CellErrorL1();
      L2Norm += SolnPtr[0][0][i].CellErrorL2();
      LMaxNorm = max(LMaxNorm,SolnPtr[0][0][i].CellErrorL1()/SolnPtr[0][0][i].CellDomain());
      Domain += SolnPtr[0][0][i].CellDomain();

      ++FinishedCell;
      Print_Progress(FinishedCell,ProgressFrequency);
    }

    /* compute final expression for norms */
    L1Norm /= Domain;
    L2Norm = sqrt(L2Norm/Domain);
    break;
  case TwoD:
    ProgressFrequency = (IP.iCell()*IP.jCell())/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=iStart(), FinishedCell=0; i<=iEnd(); ++i)
      for (int j=jStart(); j<=jEnd(); ++j){
	/* compute error in each cell */
	SolnPtr[0][j][i].ComputeReconstructionError(IP.TestFunction());
	/* update norms */
	L1Norm += SolnPtr[0][j][i].CellErrorL1();
	L2Norm += SolnPtr[0][j][i].CellErrorL2();
	LMaxNorm = max(LMaxNorm,SolnPtr[0][j][i].CellErrorL1()/SolnPtr[0][j][i].CellDomain());
	Domain += SolnPtr[0][j][i].CellDomain();

	++FinishedCell;
	Print_Progress(FinishedCell,ProgressFrequency);
      }
    /* compute final expression for norms */
    L1Norm /= Domain;
    L2Norm = sqrt(L2Norm/Domain);
    break;
  case ThreeD:
    /* Must be modified!!!! */
    ProgressFrequency = (IP.iCell()*IP.jCell()*IP.kCell())/10;
    if (ProgressFrequency == 0){
      ProgressFrequency = 10;
    }

    for (int i=iStart(), FinishedCell=0; i<=iEnd(); ++i)
      for (int j=jStart(); j<=jEnd(); ++j)
	for (int k=kStart(); k<=kEnd(); ++k){
	  /* compute error in each cell */
	  SolnPtr[k][j][i].ComputeReconstructionError(IP.TestFunction());
	  /* update norms */
	  L1Norm += SolnPtr[k][j][i].CellErrorL1();
	  L2Norm += SolnPtr[k][j][i].CellErrorL2();
	  LMaxNorm = max(LMaxNorm,SolnPtr[k][j][i].CellErrorL1()/SolnPtr[k][j][i].CellDomain());
          Domain += SolnPtr[k][j][i].CellDomain();

          ++FinishedCell;
          Print_Progress(FinishedCell,ProgressFrequency);
	}
    /* compute final expression for norms */
    L1Norm /= Domain;
    L2Norm = sqrt(L2Norm/Domain);
    break;
  }
}

/*********************************************************
* ComputeMultipleCorrelationCoefficient(void)            *
*                                                        *
* Computes the multiple correlation coefficient for each *
* computational cell based on the solution at the        *
* centroid after the reconstruction has been carried out.*
*********************************************************/
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
ComputeMultipleCorrelationCoefficient(void)
{

  vector<int> i_index;
  int CellsInOneDirection;

  CellsInOneDirection = 1 + 2*SolnPtr[0][0][0].CellRings();
  i_index.reserve(CellsInOneDirection);
  
  for (int i=iStart(); i<=iEnd(); ++i){
    /* Make Stencil */
    MakeReconstructionStencil(SolnPtr[0][0][0].CellRings(),i,i_index);
    
    /* estimate the correlation coefficient for the current cell */
    MultipleCorrelationCoefficient(*this,i_index,i);
  }
}

// PrintErrorNorms
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
void ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::
  PrintErrorNorms(void) const
{
  Print_(L1Norm);
  Print_(L2Norm);
  Print_(LMaxNorm);
}

/* Friend functions */
template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
bool operator==(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
		const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right){

  bool result = true;
  result = result && (left.N_XYZ==right.N_XYZ);
  if (result){
    switch(SpaceDimension){
    case OneD:
      for (int i=0; i<left.N_XYZ[0]; ++i)
	result = result && (left.SolnPtr[0][0][i] == right.SolnPtr[0][0][i]);
      break;
    case TwoD:
      for (int i=0; i<left.N_XYZ[0]; ++i)
	for (int j=0; j<left.N_XYZ[1]; ++j)
	   result = result && (left.SolnPtr[0][j][i] == right.SolnPtr[0][j][i]);
      break;
     case ThreeD:
       for (int i=0; i<left.N_XYZ[0]; ++i)
	 for (int j=0; j<left.N_XYZ[1]; ++j)
	   for (int k=0; k<left.N_XYZ[2]; ++k)
	     result = result && (left.SolnPtr[k][j][i] == right.SolnPtr[k][j][i]);
       
       break;
    }
  }
  return result;
}

template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
bool operator!=(const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& left,
		const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& right){
  return !(left == right);
}

template< SpaceType SpaceDimension, class GeometryType, class SolutionType> inline
std::ostream& operator<< (std::ostream& os,
			  const ComputationalDomain<SpaceDimension,GeometryType,SolutionType>& Obj){

  /* Prints out the ghost cells too */
  switch(SpaceDimension){
  case OneD:
    os << "NCx = " << Obj.N_XYZ[0] << std::endl
       << "Nghost= " << Obj.Nghost << std::endl
       << "StartIndex= " << Obj.iStart() << std::endl
       << "StopIndex= " << Obj.iEnd() << std::endl;
    for (int i=0; i<Obj.N_XYZ[0]; ++i)
      os << Obj.SolnPtr[0][0][i];
    break;
  case TwoD:
    os << "NCx = " << Obj.N_XYZ[0] << std::endl
       << "NCy = " << Obj.N_XYZ[1] << std::endl
       << "Nghost= " << Obj.Nghost << std::endl
       << "StartIndexI= " << Obj.iStart() << std::endl
       << "StopIndexI= "  << Obj.iEnd() << std::endl
       << "StartIndexJ= " << Obj.jStart() << std::endl
       << "StopIndexJ= "  << Obj.jEnd() << std::endl;
    for (int i=0; i<Obj.N_XYZ[0]; ++i)
      for (int j=0; j<Obj.N_XYZ[1]; ++j)
	os << Obj.SolnPtr[0][j][i];
    break;
  case ThreeD:
    os << "NCx = " << Obj.N_XYZ[0] << std::endl
       << "NCy = " << Obj.N_XYZ[1] << std::endl
       << "NCz = " << Obj.N_XYZ[2] << std::endl
       << "Nghost= " << Obj.Nghost << std::endl
       << "StartIndexI= " << Obj.iStart() << std::endl
       << "StopIndexI= "  << Obj.iEnd() << std::endl
       << "StartIndexJ= " << Obj.jStart() << std::endl
       << "StopIndexJ= "  << Obj.jEnd() << std::endl
       << "StartIndexK= " << Obj.kStart() << std::endl
       << "StopIndexK= "  << Obj.kEnd() << std::endl;
    for (int i=0; i<Obj.N_XYZ[0]; ++i)
      for (int j=0; j<Obj.N_XYZ[1]; ++j)
	for (int k=0; k<Obj.N_XYZ[2]; ++k)
	  os << Obj.SolnPtr[k][j][i];
    break;
  }
  return os;
}

// Static variables
template< SpaceType SpaceDimension, class GeometryType, class SolutionType>
HeaderData ComputationalDomain<SpaceDimension,GeometryType,SolutionType>::VarNames("Solution");

/*******************************************************************************
                              Set Up Function
********************************************************************************/
template<class SolutionType, class InputParameterType>
void SetUpDomain(SolutionType & SolnBlk, InputParameterType &IP)
{
  std::cout << " Set the Computational Domain:\n";
  // Set the Computation Domain
  SolnBlk.SetDomain(IP);
  
  std::cout << " Compute the intial solution:\n";
  // Compute the initial solution
  SolnBlk.SetInitialData(IP);
}

/*******************************************************************************
                              Solve Task Function
********************************************************************************/
template<class SolutionType, class InputParameterType>
void SolveTask(SolutionType & SolnBlk, InputParameterType &IP)
{

  // a switch must be introduce to distinguish between Filtering and Reconstruction
  std::cout << " Reconstruct the solution:\n";
  // Reconstruct solution
  SolnBlk.ReconstructSolution(IP);

  // Update subgrid solution
  SolnBlk.UpdateSubgridSolution();
}

/*******************************************************************************
                              Output Functions
********************************************************************************/
template<typename SolutionBlock, typename InputParameters>
void Output_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP, bool MultipleOrOnlyOneZone = true)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  if(MultipleOrOnlyOneZone){
    strcat(output_file_name, ".dat");
  } else {
    strcat(output_file_name, "_SubgridSol.dat");
  }
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension <<"D output data file(s).\n";
    cout.flush();
  }

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;

  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
  case 1:
    output_file << "st";
    break;
  case 2:
    output_file << "nd";
    break;
  case 3:
    output_file << "rd";
    break;
  default:
    output_file << "th";
  }
  output_file <<" Order " << IP.IP_SpaceDimension << "D Reconstruction Cell Solution for function "
	      << IP.Function_Type << "\"" << "\n";

  // Write solution data
  if(MultipleOrOnlyOneZone == true){
    SolnBlk.OutputSolutionCellTecplot(output_file);
  } else {
    SolnBlk.OutputSolutionCellTecplotOneZone(output_file);
  }

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();

  // ********************** Output Tecplot Piecewise Constant Solution ***************************/
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_PWC.dat");

  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension <<"D output data file(s).\n";
    cout.flush();
  }

  // Write the name of the file to the console
  cout << "\n Writing piecewise constant data to output data file: "
       << output_file_name << std::endl;

  // Write solution data
  SolnBlk.OutputPWC(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();

  // ********************** Output Smoothness Indicator Flag ***************************/

  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_SI_Flag.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open " << output_file_name << " output data file(s).\n";
    cout.flush();
  }

  // Write the name of the file to the console
  cout << "\n Writing solution smoothness indicator flag to output data file: "
       << output_file_name << std::endl;

  // Write solution data
  SolnBlk.OutputSmoothnessIndicatorFlag(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Mesh_Nodes_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_mesh.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output data file(s) for mesh geometry.\n";
    cout.flush();
  }

  cout << "\n Writing Reconstruction mesh nodes to output data file(s): ";

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;
  
  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << " Mesh nodes for grid " << IP.Grid_Type << "\"\n";

  // Write solution data
  SolnBlk.OutputMeshNodesTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Mesh_Cells_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_mesh_cells.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output data file(s) for mesh geometry.\n";
    cout.flush();
  }

  cout << "\n Writing Reconstruction mesh nodes to output data file(s): ";

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;
  
  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << " Mesh cells for grid " << IP.Grid_Type << "\"\n";

  // Write solution data
  SolnBlk.OutputMeshCellsTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Solution_Nodes_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_nodes.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output data file(s) for node solution.\n";
    cout.flush();
  }

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;
  
  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
  case 1:
    output_file << "st";
    break;
  case 2:
    output_file << "nd";
    break;
  case 3:
    output_file << "rd";
    break;
  default:
    output_file << "th";
  }
  output_file <<" Order " << IP.IP_SpaceDimension << "D Reconstruction Node Solution for function "
	      << IP.Function_Type << "\"" << "\n";

  // Write solution data
  SolnBlk.OutputSolutionNodesTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();

  // ********************** Output Smoothness Indicator Flag at Nodes ***************************/

  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_SI_Flag_nodes.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open " << output_file_name << " output data file(s).\n";
    cout.flush();
  }

  // Write the name of the file to the console
  cout << "\n Writing solution smoothness indicator flag at nodes to output data file: "
       << output_file_name << std::endl;

  // Write solution data
  SolnBlk.OutputSmoothnessIndicatorAtNodesTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();

  // ********************** Output Solution at Gauss Nodes ***************************/

  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_GQP.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open " << output_file_name << " output data file(s).\n";
    cout.flush();
  }

  // Write the name of the file to the console
  cout << "\n Writing solution at the Gauss Quadrature Points to output data file: "
       << output_file_name << std::endl;

  // Write solution data
  SolnBlk.OutputSolutionAtGaussPointsTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Full_Solution_Nodes_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_all_nodes.dat");
  
  // Open the output data file.
  
  output_file.open(output_file_name, ios::out);
  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension <<"D output data file(s) for node solution of \
                                                                     the whole domain.\n";
    cout.flush();
  }

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;  

  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
  case 1:
    output_file << "st";
    break;
  case 2:
    output_file << "nd";
    break;
  case 3:
    output_file << "rd";
    break;
  default:
    output_file << "th";
  }
  output_file <<" Order " << IP.IP_SpaceDimension << "D Reconstruction Node Solution for function "
	      << IP.Function_Type << "\"" << "\n";

  // Write solution data
  SolnBlk.OutputFullSolutionNodesTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Function_Graph_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], extension[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(extension,"_Function_");
  strcat(extension,IP.Function_Type);
  strcpy(output_file_name, prefix);
  strcat(output_file_name, extension);
  strcat(output_file_name, ".dat");
  
  // Open the output data file.
  
  output_file.open(output_file_name, ios::out);
  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output data file(s) for exact solution.\n";
    cout.flush();
  }

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;

  // Create the mesh for computation of exact solution
  SolutionBlock ExactSolnBlk;
  ExactSolnBlk = SolnBlk;	/* copy geometry */

  /* The exact solution is computed in the same subgrid points.
   For a better resolution more Subgrid points should be introduced. */
  switch(IP.IP_SpaceDimension){
  case OneD:
    for(int i=ExactSolnBlk.iStart(); i<=ExactSolnBlk.iEnd(); ++i){
      ExactSolnBlk(i).UpdateExactSubgridSolution(IP.TestFunction());
    }
    
    // Write the title.
    output_file << setprecision(14);
    output_file << "TITLE = \"" << IP.Function_Type
		<< " representation on the grid "
		<< IP.Grid_Type << "\"" << "\n";
    output_file << "VARIABLES = \"x\" \n";
    break;

  case TwoD:
    for(int i=ExactSolnBlk.iStart(); i<=ExactSolnBlk.iEnd(); ++i)
      for(int j=ExactSolnBlk.jStart(); j<=ExactSolnBlk.jEnd(); ++j){
	ExactSolnBlk(i,j).UpdateExactSubgridSolution(IP.TestFunction());
      }
    
    // Write the title.
    output_file << setprecision(14);
    output_file << "TITLE = \"" << IP.Function_Type
		<< " representation on the grid "
		<< IP.Grid_Type << "\"" << "\n";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n";
    break;

  case ThreeD:
    for(int i=ExactSolnBlk.iStart(); i<=ExactSolnBlk.iEnd(); ++i)
      for(int j=ExactSolnBlk.jStart(); j<=ExactSolnBlk.jEnd(); ++j)
	for(int k=ExactSolnBlk.kStart(); k<=ExactSolnBlk.kEnd(); ++k){
	  ExactSolnBlk(i,j,k).UpdateExactSubgridSolution(IP.TestFunction());
	}
    
    // Write the title.
    output_file << setprecision(14);
    output_file << "TITLE = \"" << IP.Function_Type
		<< " representation on the grid "
		<< IP.Grid_Type << "\"" << "\n";
    output_file << "VARIABLES = \"x\" \\ \n"
		<< "\"y\" \\ \n"
		<< "\"z\" \\ \n";
    break;
  }

  SolutionBlock::VarNames.PrintHeaderTecplot(output_file);
  /* output dataset auxiliary data */
  TecplotDatasetAuxData::PrintAuxData(output_file);


  // Write solution data
  ExactSolnBlk.OutputSolutionCellTecplotOneZone(output_file,false);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

// Write the error norms (L1, L2 and LMax) to file
template<typename SolutionBlock, typename InputParameters>
void Output_Error_Norms_Tecplot(SolutionBlock & SolnBlk, InputParameters &IP, bool & Title_Error_Norms)
{
  
  /* open the stream for outputing */
  int i;
  char prefix[256], extension[256], output_file_name[256];
  char output_directory[30];
  char order[3];
  ofstream output_file;

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  switch(IP.Reconstruction_Order){
  case 1:
    strcpy(order,"1");
    break;
  case 2:
    strcpy(order,"2");
    break;
  case 3:
    strcpy(order,"3");
    break;
  case 4:
    strcpy(order,"4");
    break;
  case 5:
    strcpy(order,"5");
    break;
  case 6:
    strcpy(order,"6");
    break;
  case 7:
    strcpy(order,"7");
    break;
  default:
    strcpy(order,"0");
  }

  strcpy(extension,"_Function_");
  strcat(extension,IP.Function_Type);
  strcat(extension,"_LNorms_Order");
  strcat(extension,order);

  // Set output_file_name;
  strcpy(output_file_name, prefix);
  strcat(output_file_name, extension);
  strcat(output_file_name, ".dat");

  // Open the output data files.

  if (Title_Error_Norms){

    // L Norms
    output_file.open(output_file_name, ios::trunc);
    if (output_file.bad()) {
      cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	   << " Unable to open Reconstruct" << IP.IP_SpaceDimension <<"D output data file(s) for error norms.\n";
      cout.flush();
    }

    // Write the name of the file to the console
    std::cout << "\n: "<<  output_file_name;  

    // write the header
    switch(IP.IP_SpaceDimension){
    case OneD:
      /* 1D */
      output_file << "TITLE = \"" << " LNorms for function "
		  << IP.Function_Type << " reconstructed on the grid "
		  << IP.Grid_Type <<"\""<< "\n"
		  << "#    #Nodes =  pow(TotalNumberOfNodes,1/" << IP.IP_SpaceDimension << ")\n" 
		  << "VARIABLES = \"N<sub>Nodes</sub>\" \n " 
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>1</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>2</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><math>%</math></sub> \" \n";

      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE \n";
      break;

    case TwoD:
      /* 2D */
      output_file << "TITLE = \"" << " LNorms for function "
		  << IP.Function_Type << " reconstructed on the grid "
		  << IP.Grid_Type <<"\""<< "\n"
		  << "#    #Nodes =  pow(TotalNumberOfNodes,1/" << IP.IP_SpaceDimension << ")\n" 
		  << "VARIABLES = \"1/N<sup>1/2</sup><sub>Nodes</sub>\" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>1</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>2</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><math>%</math></sub> \" \n";

      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE \n";
      break;

    case ThreeD:
      /* 3D */
      output_file << "TITLE = \"" << " LNorms for function "
		  << IP.Function_Type << " reconstructed on the grid "
		  << IP.Grid_Type <<"\""<< "\n"
		  << "#    #Nodes =  pow(TotalNumberOfNodes,1/" << IP.IP_SpaceDimension << ")\n" 
		  << "VARIABLES = \"1/N<sup>1/3</sup><sub>Nodes</sub>\" \n " 
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>1</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><greek>2</greek></sub> \" \n "
		  << " \"<math>o</math><greek>E</greek><math>o</math><sub><math>%</math></sub> \" \n";

      /* output dataset auxiliary data */
      TecplotDatasetAuxData::PrintAuxData(output_file);

      output_file << "ZONE \n";
      break;
    }

    // Make sure that the title is written only once
    Title_Error_Norms = false;
  }
  else {
    output_file.open(output_file_name, ios::app);
    if (output_file.bad()) {
      cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	   << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output  error norms data file(s) \
                                                                        for adding data.\n";
      cout.flush();
    }
  }
  
  // Write Norms
  output_file << setprecision(14);
  output_file << " " << pow(1.0*SolnBlk.NumberOfNodes(),1.0/IP.IP_SpaceDimension) << "\t"
    /* Obs. The mesh size "DeltaX" is proportional with 1.0/pow(#Nodes,1/IP.IP_SpaceDimension) */
	      << SolnBlk.CellL1Norm()   << "\t" 
	      << SolnBlk.CellL2Norm()   << "\t" 
	      << SolnBlk.CellLMaxNorm() << "\n";
  output_file.unsetf(ios::scientific);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

template<typename SolutionBlock, typename InputParameters>
void Output_Tecplot_Stencil_Reconstruction(SolutionBlock & SolnBlk, InputParameters &IP)
{
  /* open the stream for outputing */

  int i;
  char prefix[256], output_file_name[256], output_directory[30];
  ofstream output_file;    

  // Determine prefix of output data file names.
  
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name) ) break;
  } // endwhile
  prefix[i] = '\0';
  
  // Determine output data file name.
  strcpy(output_file_name, prefix);
  strcat(output_file_name, "_StencilRec.dat");
  
  // Open the output data file.
  output_file.open(output_file_name, ios::out);

  if (output_file.bad()) {
    cout << "\n Reconstruct" << IP.IP_SpaceDimension << "D ERROR:" 
	 << " Unable to open Reconstruct" << IP.IP_SpaceDimension << "D output data file(s) for stencil subgrid solution.\n";
    cout.flush();
  }

  // Write the name of the file to the console
  std::cout << output_file_name << std::endl;
  
  // Write the title.
  output_file << setprecision(14);
  output_file << "TITLE = \"" << IP.Reconstruction_Order;
  switch (IP.Reconstruction_Order){
  case 1:
    output_file << "st";
    break;
  case 2:
    output_file << "nd";
    break;
  case 3:
    output_file << "rd";
    break;
  default:
    output_file << "th";
  }
  output_file <<" Order " << IP.IP_SpaceDimension << "D Reconstruction Solution in every stencil for function "
	      << IP.Function_Type << "\"" << "\n";

  // Write solution data
  SolnBlk.OutputReconstructedSolutionTecplot(output_file);

  // Close the output data file.
  output_file << setprecision(6);  
  output_file.close();
}

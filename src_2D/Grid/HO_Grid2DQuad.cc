/*!\file HO_Grid2DQuad.cc
   \brief Source file initializing/implementing member variables/functions that belong to classes defined in HO_Grid2DQuad.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuad.h"	// Include 2D high-order block grid

// ===== Member variables =====

/*!
 * This switch is used to determine whether the geometric properties
 * (i.e. centroid, area etc.) of the cells closed to boundaries are
 * determined using a high-order or a low-order representation of the boundary.
 * If ON, the boundary splines are used to compute the geometric properties (i.e. high-order accuracy).
 * If OFF, a linear representation between 2 consecutive nodes is considered (i.e. 2nd-order accuracy).
 */
int Grid2D_Quad_Block_HO::HighOrderBoundaryRepresentation = OFF;

/*!
 * This switch is used to determine whether the smoothing 
 * subroutine is applied or not. 
 * This flag acts globally, that is, on all quad blocks.
 * If ON, the smoother is applied.
 * If OFF, the smoother is not applied so a call to the smoothing subroutine leaves the mesh unchanged.
 */
int Grid2D_Quad_Block_HO::Smooth_Quad_Block_Flag = ON;

// ===== Member functions =====

/*!
 * Default constructor.
 */
Grid2D_Quad_Block_HO::Grid2D_Quad_Block_HO(void)
  : Integration(this),
    NNi(0), INl(0), INu(0), NNj(0), JNl(0), JNu(0),
    NCi(0), ICl(0), ICu(0), NCj(0), JCl(0), JCu(0),
    Nghost(0),
    Node(NULL), Cell(NULL),
    BCtypeN(NULL), BCtypeS(NULL), BCtypeE(NULL), BCtypeW(NULL),
    BndNorthSpline(), BndSouthSpline(), BndEastSpline(), BndWestSpline(),
    BndNorthSplineInfo(NULL), BndSouthSplineInfo(NULL),
    BndEastSplineInfo(NULL), BndWestSplineInfo(NULL),
    SminN(ZERO), SmaxN(ZERO), SminS(ZERO), SmaxS(ZERO), 
    SminE(ZERO), SmaxE(ZERO), SminW(ZERO), SmaxW(ZERO),
    StretchI(0), StretchJ(0), BetaI(ONE), TauI(ONE),
    BetaJ(ONE), TauJ(ONE),
    OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1),
    // Initialize mesh update flags to OFF (i.e. no update scheduled)
    InteriorMeshUpdate(OFF), GhostCellsUpdate(OFF), CornerGhostCellsUpdate(OFF)
{
  // 
}

/*!
 * Copy constructor. It is declared private
 */
Grid2D_Quad_Block_HO::Grid2D_Quad_Block_HO(const Grid2D_Quad_Block_HO &G)
  :Integration(this),
   NNi(0), INl(0), INu(0), NNj(0), JNl(0), JNu(0),
   NCi(0), ICl(0), ICu(0), NCj(0), JCl(0), JCu(0),
   Nghost(0),
   Node(NULL), Cell(NULL),
   BCtypeN(NULL), BCtypeS(NULL), BCtypeE(NULL), BCtypeW(NULL),
   BndNorthSpline(), BndSouthSpline(), BndEastSpline(), BndWestSpline(),
   BndNorthSplineInfo(NULL), BndSouthSplineInfo(NULL),
   BndEastSplineInfo(NULL), BndWestSplineInfo(NULL),
   SminN(ZERO), SmaxN(ZERO), SminS(ZERO), SmaxS(ZERO), 
   SminE(ZERO), SmaxE(ZERO), SminW(ZERO), SmaxW(ZERO),
   StretchI(0), StretchJ(0),
   BetaI(ONE), TauI(ONE),
   BetaJ(ONE), TauJ(ONE),
   OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1),
   // Initialize mesh update flags to OFF (i.e. no update scheduled)
   InteriorMeshUpdate(OFF), GhostCellsUpdate(OFF), CornerGhostCellsUpdate(OFF)
{
  int Ni, Nj;
  int i,j;

  // allocate memory for the new container
  Ni = G.NCi - 2*G.Nghost;
  Nj = G.NCj - 2*G.Nghost;
  allocate(Ni,Nj,G.Nghost);

  // Set the grid values by copying from the grid block G.
  if (G.Node != NULL) {

    // Copy the node locations of the grid block G.
    for (j  = G.JNl-G.Nghost; j <= G.JNu+G.Nghost ; ++j ) {
      for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost ; ++i ) {
	Node[i][j].X = G.Node[i][j].X;
      } /* endfor */
    } /* endfor */
    
    // Copy the cell values of grid block G.
    for ( j = G.JCl-G.Nghost; j <= G.JCu+G.Nghost ; ++j) {
      for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
	Cell[i][j] = G.Cell[i][j];
      } /* endfor */
    } /* endfor */

    // Copy the boundary condition type info of grid block G.
    for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
      BCtypeN[i] = G.BCtypeN[i];
      BCtypeS[i] = G.BCtypeS[i];
    } /* endfor */
    for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
      BCtypeE[j] = G.BCtypeE[j];
      BCtypeW[j] = G.BCtypeW[j];
    } /* endfor */
  } /* endif */


  // Copy boundary spline info of grid block G.
  if (G.BndNorthSpline.np != 0) {
    BndNorthSpline = G.BndNorthSpline;
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
  } /* endif */

  if (G.BndSouthSpline.np != 0) {
    BndSouthSpline = G.BndSouthSpline;
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
  } /* endif */
  
  if (G.BndEastSpline.np != 0) {
    BndEastSpline = G.BndEastSpline;
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
  } /* endif */
  
  if (G.BndWestSpline.np != 0) {
    BndWestSpline = G.BndWestSpline;
  } else if (BndWestSpline.np != 0) {
    BndWestSpline.deallocate();
  } /* endif */

  // Copy boundary spline pathlength info of grid block G.
  SminN = G.SminN;
  SmaxN = G.SmaxN;
  SminS = G.SminS;
  SmaxS = G.SmaxS;
  SminE = G.SminE;
  SmaxE = G.SmaxE;
  SminW = G.SminW;
  SmaxW = G.SmaxW;
  
  // Copy node stretching info of grid block G.
  StretchI = G.StretchI;
  BetaI = G.BetaI;
  TauI = G.TauI;
  StretchJ = G.StretchJ;
  BetaJ = G.BetaJ;
  TauJ = G.TauJ;
  OrthogonalN = G.OrthogonalN;
  OrthogonalS = G.OrthogonalS;
  OrthogonalE = G.OrthogonalE;
  OrthogonalW = G.OrthogonalW;

  /* No mesh update is required for this operation
     since all geometric values are copied directly. */
  Reset_Mesh_Update_Flags();
}

/*!
 * Allocate memory.
 *
 * \param Ni number of cells in i-direction
 * \param Nj number of cells in j-direction
 * \param Ng number of ghost cells
 * \param HighestRecOrder the highest reconstruction order used during the computation.
 */
void Grid2D_Quad_Block_HO::allocate(const int &Ni, const int &Nj, const int &Ng, const int &HighestRecOrder) {
  int i,j;

  // Check conditions
  assert( Ni > 1 && Nj > 1 && Ng >= 2 && HighestRecOrder >= 0);

  // Check if the new required memory has dimensions different than the currently allocated ones
  if ( Ni != (NCi-2*Nghost) || Nj != (NCj-2*Nghost) || Ng != Nghost ){

    // free the memory if there is memory allocated
    deallocate();
    
    // allocate new memory 
    NNi = Ni+2*Ng+1; INl = Ng; INu = Ni+Ng;
    NNj = Nj+2*Ng+1; JNl = Ng; JNu = Nj+Ng;
    NCi = Ni+2*Ng;   ICl = Ng; ICu = Ni+Ng-1;
    NCj = Nj+2*Ng;   JCl = Ng; JCu = Nj+Ng-1;
    Nghost = Ng;
    
    Node = new Node2D_HO*[NNi];
    for ( i = 0; i <= NNi-1 ; ++i ) Node[i] = new Node2D_HO[NNj];
    Cell = new Cell2D_HO*[NCi];
    for ( i = 0; i <= NCi-1 ; ++i ) Cell[i] = new Cell2D_HO[NCj];

    // allocate memory for the container of geometric moments
    for ( i = 0; i <NCi ; ++i ){
      for ( j = 0; j <NCj ; ++j ){
	Cell[i][j].SetGeomCoeffContainerSize(HighestRecOrder);       
      }
    }

    BCtypeN = new int[NCi]; BCtypeS = new int[NCi];
    BCtypeE = new int[NCj]; BCtypeW = new int[NCj];

    // Complete memory allocation
    return;
    
  } else if (HighestRecOrder != Cell[0][0].GeomCoeff().RecOrder()){
    // Check if the highest reconstruction order is different than the current one.
    
    // re-allocate memory for the container of geometric moments
    for ( i = 0; i <NCi ; ++i ){
      for ( j = 0; j <NCj ; ++j ){
	Cell[i][j].SetGeomCoeffContainerSize(HighestRecOrder);       
      }
    }    
  }
}

/*
 * Deallocate memory.
 */
void Grid2D_Quad_Block_HO::deallocate(void) {
  int i;

  // Deallocate nodes
  if (Node != NULL){
    for ( i = 0; i <= NNi-1 ; ++i ) {
      delete []Node[i]; Node[i] = NULL;
    } /* endfor */
    delete []Node; Node = NULL;
  }

  // Deallocate cells
  if (Cell != NULL){
    for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []Cell[i]; Cell[i] = NULL;
    } /* endfor */
    delete []Cell; Cell = NULL; 
  }

  // Deallocate North boundary conditions
  if (BCtypeN != NULL){
    delete []BCtypeN; BCtypeN = NULL;
  }
  // Deallocate South boundary conditions
  if (BCtypeS != NULL){
    delete []BCtypeS; BCtypeS = NULL; 
  }
  // Deallocate East boundary conditions
  if (BCtypeE != NULL){
    delete []BCtypeE; BCtypeE = NULL;
  }
  // Deallocate West boundary conditions
  if (BCtypeW != NULL){
    delete []BCtypeW; BCtypeW = NULL;
  }

  // Deallocate boundary splines
  BndNorthSpline.deallocate(); BndSouthSpline.deallocate();
  BndEastSpline.deallocate(); BndWestSpline.deallocate();

  // Deallocate boundary spline info
  deallocateBndSplineInfo();

  // Reset mesh indexes
  NNi = 0; INl = 0; INu = 0; NNj = 0; JNl = 0; JNu = 0;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0;
  Nghost = 0;
  StretchI = 0; StretchJ = 0; BetaI = ONE; TauI = ONE;
  BetaJ = ONE; TauJ = ONE;
  OrthogonalN = 1; OrthogonalS = 1; OrthogonalE = 1; OrthogonalW = 1;
}

/*!
 * Assignment operator =
 */
Grid2D_Quad_Block_HO& Grid2D_Quad_Block_HO::operator=(const Grid2D_Quad_Block_HO &Grid) {

  // !!! If the LHS grid block already has objects assigned, these are going to be deleted.
  // Handle self-assignment:
  if (this == & Grid) return *this;

  int Ni, Nj;
  int i,j;

  // re-allocate memory if there isn't enough
  Ni = Grid.NCi - 2*Grid.Nghost;
  Nj = Grid.NCj - 2*Grid.Nghost;
  allocate(Ni,Nj,Grid.Nghost);

  // Set the grid values by copying from the grid block Grid.
  if (Grid.Node != NULL) {

    // Copy the node locations of the grid block Grid.
    for (j  = Grid.JNl-Grid.Nghost; j <= Grid.JNu+Grid.Nghost ; ++j ) {
      for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost ; ++i ) {
	Node[i][j].X = Grid.Node[i][j].X;
      } /* endfor */
    } /* endfor */
    
    // Copy the cell values of grid block Grid.
    for ( j = Grid.JCl-Grid.Nghost; j <= Grid.JCu+Grid.Nghost ; ++j) {
      for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
	Cell[i][j] = Grid.Cell[i][j];
      } /* endfor */
    } /* endfor */

    // Copy the boundary condition type info of grid block Grid.
    for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
      BCtypeN[i] = Grid.BCtypeN[i];
      BCtypeS[i] = Grid.BCtypeS[i];
    } /* endfor */
    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
      BCtypeE[j] = Grid.BCtypeE[j];
      BCtypeW[j] = Grid.BCtypeW[j];
    } /* endfor */
  } /* endif */


  // Copy boundary spline info of grid block Grid.
  if (Grid.BndNorthSpline.np != 0) {
    BndNorthSpline = Grid.BndNorthSpline;
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
  } /* endif */

  if (Grid.BndSouthSpline.np != 0) {
    BndSouthSpline = Grid.BndSouthSpline;
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
  } /* endif */
  
  if (Grid.BndEastSpline.np != 0) {
    BndEastSpline = Grid.BndEastSpline;
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
  } /* endif */
  
  if (Grid.BndWestSpline.np != 0) {
    BndWestSpline = Grid.BndWestSpline;
  } else if (BndWestSpline.np != 0) {
    BndWestSpline.deallocate();
  } /* endif */

  // Copy boundary spline pathlength info of grid block Grid.
  SminN = Grid.SminN;
  SmaxN = Grid.SmaxN;
  SminS = Grid.SminS;
  SmaxS = Grid.SmaxS;
  SminE = Grid.SminE;
  SmaxE = Grid.SmaxE;
  SminW = Grid.SminW;
  SmaxW = Grid.SmaxW;
  
  // Copy node stretching info of grid block Grid.
  StretchI = Grid.StretchI;
  BetaI = Grid.BetaI;
  TauI = Grid.TauI;
  StretchJ = Grid.StretchJ;
  BetaJ = Grid.BetaJ;
  TauJ = Grid.TauJ;
  OrthogonalN = Grid.OrthogonalN;
  OrthogonalS = Grid.OrthogonalS;
  OrthogonalE = Grid.OrthogonalE;
  OrthogonalW = Grid.OrthogonalW;

  /* No mesh update is required for this operation
     since all geometric values are copied directly. */
  Reset_Mesh_Update_Flags();
}

/*!
 * Calculate influence coefficients 
 * based on bilinear interpolation.
 */
double Grid2D_Quad_Block_HO::dNdC(const int ii, const int jj, 
				  const int cell_orientation) const {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
    eta1, zeta1, eta2, zeta2, eta, zeta;
  double interpolation_coefficient;
  x  = Node[ii][jj].X.x;
  y  = Node[ii][jj].X.y;
  ax = Cell[ii-1][jj-1].Xc.x;
  bx = Cell[ii-1][jj  ].Xc.x - Cell[ii-1][jj-1].Xc.x; 
  cx = Cell[ii  ][jj-1].Xc.x - Cell[ii-1][jj-1].Xc.x; 
  dx = Cell[ii  ][jj  ].Xc.x + Cell[ii-1][jj-1].Xc.x -
    Cell[ii-1][jj  ].Xc.x - Cell[ii  ][jj-1].Xc.x;
  ay = Cell[ii-1][jj-1].Xc.y;
  by = Cell[ii-1][jj  ].Xc.y - Cell[ii-1][jj-1].Xc.y; 
  cy = Cell[ii  ][jj-1].Xc.y - Cell[ii-1][jj-1].Xc.y; 
  dy = Cell[ii  ][jj  ].Xc.y + Cell[ii-1][jj-1].Xc.y -
    Cell[ii-1][jj  ].Xc.y - Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by+dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) {
      zeta1 = -cc/bb;
    } else { 
      zeta1 = -cc/sgn(bb)*(TOLER*TOLER); 
    } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } else { 
      eta1 = HALF;
    } 
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { zeta1 = -HALF*bb/aa;
    } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    }
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) {
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else {
      zeta2 = -HALF*bb/aa;
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) { 
      eta2 = -ONE;
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  // Determine the coefficient or weight for cell orientation.
  switch(cell_orientation) {
  case NORTH_EAST : 
    interpolation_coefficient = zeta*eta;
    break;
  case NORTH_WEST : 
    interpolation_coefficient = zeta - zeta*eta;
    break;
  case SOUTH_EAST : 
    interpolation_coefficient = eta - zeta*eta;
    break;
  case SOUTH_WEST : 
    interpolation_coefficient = ONE - zeta - eta + zeta*eta;
    break;
  default:
    cout <<"\n Improper Orientation in Grid2D_Quad::dNdC!\n"; 
    interpolation_coefficient = -ONE;
    break;
  }
  return interpolation_coefficient;
}

/*!
 * Calculate influence coefficients 
 * based on bilinear interpolation. 
 */
//  dFacedC:
//  
//  We know that we can write each face variable, uf, as a geometric
//  weight of the cell-centred value, uc: uf = intp uc
//  
//  This function determine what this "intp" is.
//  
//  For diamond-path integration: (which is how we find the gradients of
//  the face variables but we are sloppy and reuse the term "diamond" also
//  to apply to how we find the variables themselves) to find the face
//  variables we first find the nodal values using an interpolation and
//  then use those nodal values together with the cell-centred variables
//  in a second (third, really, after two nodes) interpolation to find the
//  face variable.
//  
//  For example, suppose we wanted to find density at the north face of
//  cell (i, j). We would first find the value at node NE using the
//  values from the four neighbouring cells of node NE. Then we would find
//  the value at node NW. Then we would apply a third interpolation using
//  the values at NW, NE, the centre of (i, j) and the centre of (i, j+1).
//  
//  In a simplifying assumption, we ignore the third interpolation and
//  simply assume that the face value is the average of the node values. 
//  
//  In what follows, then, "left_node_weight" and "right_node_weight" 
//  are the weights used to obtain the nodal values. If the cell only 
//  contributes to one of the nodal values then one of "left_node_weight" 
//  or "right_node_weight" is zero. "intp" then is simply the average of 
//  "left_node_weight" and "right_node_weight".
//  
//  The "source" cell is in the direction "cell_orientation" from 
//  cell (i, j).
//  
//  For example, suppose that 
//
//  face = North 
//  cell_orientation = East
//
//  Now we wish to know how the cell to the east of (i, j) influences
//  the value on the north face of (i, j). In this case, the cell to
//  the east influences the value at the node to the right of the face
//  value but does not influence the value at the node to the left of the
//  face.  So "left_node_weight" is zero. Also the node to the right of 
//  the face is the NW node of the cell to the east.
//
void Grid2D_Quad_Block_HO::dFacedC(double &left_node_weight, 
				   double &right_node_weight, 
				   const int i, const int j, const int face, 
				   const int cell_orientation) const {
  int lnodei = 0, lnodej = 0; //  left node (not cell) index
  int rnodei = 0, rnodej = 0; // right node (not cell) index

  switch (face){
  case NORTH:
    lnodei = i  ; lnodej = j+1;
    rnodei = i+1; rnodej = j+1;
    if (cell_orientation == CENTER) {    
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);
    } else if (cell_orientation == NORTH) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);
    } else if (cell_orientation == EAST) {
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);
    } else if (cell_orientation == WEST){ 
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);      
      right_node_weight = ZERO;                   
    } else if (cell_orientation == NORTH_WEST) {
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);      
      right_node_weight = ZERO;   
    } else if (cell_orientation == NORTH_EAST) {
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);    
    } /* endif */
    break;
      
  case EAST:
    lnodei = i+1; lnodej = j+1;
    rnodei = i+1; rnodej = j  ;
    if (cell_orientation == CENTER) {  
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST); 
    } else if (cell_orientation == EAST) {    
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST); 
    } else if (cell_orientation == NORTH) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH) {           
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);       
    } else if (cell_orientation == NORTH_EAST) {     
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH_EAST) {           
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);       
    } /* endif */
    break;
  
  case SOUTH:
    lnodei = i+1; lnodej = j;
    rnodei = i  ; rnodej = j;
    if (cell_orientation == CENTER){ 
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);
    } else if (cell_orientation == SOUTH) {     
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);
    } else if (cell_orientation == EAST) { 
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);      
      right_node_weight = ZERO;
    } else if (cell_orientation == WEST){    
      left_node_weight = ZERO; 
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);                              
    } else if (cell_orientation == SOUTH_EAST) { 
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);      
      right_node_weight = ZERO;
    } else if (cell_orientation == SOUTH_WEST) {    
      left_node_weight = ZERO; 
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);      
    } /* endif */
    break;
   
  case WEST:
    lnodei = i; lnodej = j  ;
    rnodei = i; rnodej = j+1;
    if (cell_orientation == CENTER) {  
      left_node_weight = dNdC(lnodei, lnodej, NORTH_EAST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_EAST);  
    } else if (cell_orientation == WEST) {  
      left_node_weight = dNdC(lnodei, lnodej, NORTH_WEST);
      right_node_weight = dNdC(rnodei, rnodej, SOUTH_WEST);  
    } else if (cell_orientation == SOUTH) {           
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_EAST);
      right_node_weight = ZERO;
    } else if (cell_orientation == NORTH) {   
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_EAST);
    } else if (cell_orientation == SOUTH_WEST) {           
      left_node_weight = dNdC(lnodei, lnodej, SOUTH_WEST);
      right_node_weight = ZERO;
    } else if (cell_orientation == NORTH_WEST) {   
      left_node_weight = ZERO;
      right_node_weight = dNdC(rnodei, rnodej, NORTH_WEST);
    } /* endif */
    break;     
  }
  if (left_node_weight < ZERO || right_node_weight < ZERO) {
    cout <<"\n Invalid node weightings in Grid2D_Quad::dFacedC!\n"; 
  }
}

/*!
 * Influence coefficients for diamond path gradient
 * reconstruction of face solution gradients.
 */
//  dDiamondPathdC:
//  
//  Finds the weights used to determine the gradient of a variable on the
//  cell face from a cell-centred value. That is, if ux is the gradient of
//  velocity at a cell face and u is the velocity at one cell centre then
//  we approximate ux by:
//  
//    ux = ( d_dWdx_dW )( u ) + contributions from other cells
//  
//  This function finds d_dWdx_dW and the corresponding d_dWdy_dW. 
//
void Grid2D_Quad_Block_HO::dDiamondPathdC(double &d_dWdx_dW, double &d_dWdy_dW, 
					  const int i, const int j, const int face, 
					  const double &face_left_node_weight, 
					  const double &face_right_node_weight, 
					  const int cell_orientation) const {
  double AREA = 0.0;
  Vector2D norm[4];
  double  dWnNWdWc = 0.0, dWnNEdWc = 0.0,  dWnSWdWc = 0.0, dWnSEdWc = 0.0;
 
  switch(face){
    /*************** NORTH ****************************/
  case NORTH: 
    dWnNWdWc = face_left_node_weight;
    dWnNEdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x = nodeNE(i,j).X.y - Cell[i][j].Xc.y;
    norm[0].y = -( nodeNE(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = Cell[i][j+1].Xc.y - nodeNE(i,j).X.y;
    norm[1].y = -(Cell[i][j+1].Xc.x - nodeNE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  nodeNW(i,j).X.y - Cell[i][j+1].Xc.y;
    norm[2].y = -(nodeNW(i,j).X.x - Cell[i][j+1].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = Cell[i][j].Xc.y - nodeNW(i,j).X.y ;
    norm[3].y = -( Cell[i][j].Xc.x - nodeNW(i,j).X.x);
    AREA =  HALF*(fabs((nodeNE(i,j).X-Cell[i][j].Xc)^
		       (nodeNW(i,j).X-Cell[i][j].Xc)) +
		  fabs((nodeNW(i,j).X-Cell[i][j+1].Xc)^
		       (nodeNE(i,j).X-Cell[i][j+1].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnNEdWc)*norm[0].x+dWnNEdWc* norm[1].x+ 
			dWnNWdWc* norm[2].x+ (ONE+dWnNWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnNEdWc)*norm[0].y+dWnNEdWc*norm[1].y+ 
			dWnNWdWc* norm[2].y+ (ONE+dWnNWdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == NORTH) {
      d_dWdx_dW = HALF*(dWnNEdWc*norm[0].x + (ONE+dWnNEdWc)* norm[1].x + 
			(ONE+dWnNWdWc)* norm[2].x + dWnNWdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*norm[0].y + (ONE+dWnNEdWc)* norm[1].y + 
			(ONE+dWnNWdWc)* norm[2].y + dWnNWdWc*norm[3].y)/AREA;  
    } else if (cell_orientation == EAST || cell_orientation == NORTH_EAST) {
      d_dWdx_dW = HALF*(dWnNEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if (cell_orientation == WEST || cell_orientation == NORTH_WEST) {
      d_dWdx_dW = HALF*(dWnNWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnNWdWc*(norm[2].y+norm[3].y))/AREA;  
    } /* endif */
    break;
    
    /*************** EAST ****************************/
  case EAST:
    dWnNEdWc = face_left_node_weight;
    dWnSEdWc = face_right_node_weight; 
    //  normal vector of the SE side of a diamond 
    norm[0].x =  Cell[i+1][j].Xc.y - nodeSE(i,j).X.y;
    norm[0].y = -(Cell[i+1][j].Xc.x - nodeSE(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = nodeNE(i,j).X.y -  Cell[i+1][j].Xc.y ;
    norm[1].y = -(nodeNE(i,j).X.x -  Cell[i+1][j].Xc.x );
    //  normal vector of the NW side of a diamond 
    norm[2].x =   Cell[i][j].Xc.y - nodeNE(i,j).X.y ;
    norm[2].y = -(Cell[i][j].Xc.x - nodeNE(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = nodeSE(i,j).X.y - Cell[i][j].Xc.y;
    norm[3].y = -(nodeSE(i,j).X.x - Cell[i][j].Xc.x);
    AREA =  HALF*(fabs((nodeNE(i,j).X-Cell[i+1][j].Xc)^
		       (nodeSE(i,j).X-Cell[i+1][j].Xc)) +
		  fabs((nodeSE(i,j).X-Cell[i][j].Xc)^
		       (nodeNE(i,j).X-Cell[i][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc*norm[0].x+ dWnNEdWc*norm[1].x+ 
			(ONE+ dWnNEdWc)*norm[2].x+ (ONE+dWnSEdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*norm[0].y+ dWnNEdWc* norm[1].y+ 
			(ONE+ dWnNEdWc)*norm[2].y+ (ONE+dWnSEdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == EAST) {
      d_dWdx_dW = HALF*((ONE+dWnSEdWc)*norm[0].x + (ONE+dWnNEdWc)*norm[1].x + 
			dWnNEdWc* norm[2].x + dWnSEdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSEdWc)* norm[0].y + (ONE+dWnNEdWc)*norm[1].y + 
			dWnNEdWc*norm[2].y + dWnSEdWc*norm[3].y)/AREA;  
    } else if (cell_orientation == NORTH || cell_orientation == NORTH_EAST) {
      d_dWdx_dW = HALF*(dWnNEdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if (cell_orientation == SOUTH || cell_orientation == SOUTH_EAST) {
      d_dWdx_dW = HALF*(dWnSEdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*(norm[0].y+norm[3].y))/AREA;  
    } /* endif */
    break;

    /*************** SOUTH ****************************/
  case SOUTH:
    dWnSEdWc = face_left_node_weight;
    dWnSWdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x =  nodeSE(i,j).X.y - Cell[i][j-1].Xc.y;
    norm[0].y = -(nodeSE(i,j).X.x - Cell[i][j-1].Xc.x  );
    //  normal vector of the NE side of a diamond 
    norm[1].x = Cell[i][j].Xc.y -  nodeSE(i,j).X.y;
    norm[1].y = -(Cell[i][j].Xc.x -  nodeSE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =   nodeSW(i,j).X.y - Cell[i][j].Xc.y ;
    norm[2].y = -(nodeSW(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  Cell[i][j-1].Xc.y - nodeSW(i,j).X.y;
    norm[3].y = -(Cell[i][j-1].Xc.x- nodeSW(i,j).X.x);
    AREA =  HALF*(fabs((nodeSE(i,j).X-Cell[i][j-1].Xc)^
		       (nodeSW(i,j).X-Cell[i][j-1].Xc)) +
		  fabs((nodeSE(i,j).X-Cell[i][j].Xc)^
		       (nodeSW(i,j).X-Cell[i][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc*norm[0].x+ (ONE+dWnSEdWc)*norm[1].x+ 
			(ONE+ dWnSWdWc)*norm[2].x+ (dWnSWdWc)*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*norm[0].y+ (ONE+dWnSEdWc)*norm[1].y+ 
			(ONE+ dWnSWdWc)*norm[2].y+ (dWnSWdWc)*norm[3].y)/AREA;  
    } else if( cell_orientation == SOUTH) {
      d_dWdx_dW = HALF*((ONE+dWnSEdWc)*norm[0].x + dWnSEdWc*norm[1].x + 
			dWnSWdWc*norm[2].x + (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSEdWc)*norm[0].y + dWnSEdWc*norm[1].y + 
			dWnSWdWc*norm[2].y + (ONE+dWnSWdWc)* norm[3].y)/AREA;  
    } else if (cell_orientation == EAST || cell_orientation == SOUTH_EAST) {
      d_dWdx_dW = HALF*(dWnSEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if (cell_orientation == WEST || cell_orientation == SOUTH_WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*(norm[2].y+norm[3].y))/AREA;  
    } /* endif */
    break;
    
    /*************** WEST ****************************/
  case WEST:
    dWnSWdWc = face_left_node_weight;
    dWnNWdWc = face_right_node_weight;
    //  normal vector of the SE side of a diamond 
    norm[0].x =   Cell[i][j].Xc.y - nodeSW(i,j).X.y;
    norm[0].y = -(Cell[i][j].Xc.x - nodeSW(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x =  nodeNW(i,j).X.y - Cell[i][j].Xc.y;
    norm[1].y = -(nodeNW(i,j).X.x - Cell[i][j].Xc.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  Cell[i-1][j].Xc.y - nodeNW(i,j).X.y ;
    norm[2].y = -( Cell[i-1][j].Xc.x - nodeNW(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  nodeSW(i,j).X.y - Cell[i-1][j].Xc.y;
    norm[3].y = -(nodeSW(i,j).X.x - Cell[i-1][j].Xc.x);
    AREA =  HALF*(fabs((nodeNW(i,j).X-Cell[i][j].Xc)^
		       (nodeSW(i,j).X-Cell[i][j].Xc)) +
		  fabs((nodeNW(i,j).X-Cell[i-1][j].Xc)^
		       (nodeSW(i,j).X-Cell[i-1][j].Xc)));
    if(cell_orientation == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnSWdWc)*norm[0].x+ (ONE+dWnNWdWc)*norm[1].x+ 
			dWnNWdWc* norm[2].x+ dWnSWdWc* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSWdWc)*norm[0].y+ (ONE+dWnNWdWc)*norm[1].y+ 
			dWnNWdWc* norm[2].y+ dWnSWdWc* norm[3].y)/AREA;  
    } else if (cell_orientation == WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*norm[0].x + dWnNWdWc*norm[1].x + 
			(ONE+dWnNWdWc)*norm[2].x+ (ONE+dWnSWdWc)*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*norm[0].y + dWnNWdWc*norm[1].y + 
			(ONE+dWnNWdWc)* norm[2].y+ (ONE+dWnSWdWc)*norm[3].y)/AREA;    
    } else if (cell_orientation == NORTH || cell_orientation == NORTH_WEST) {
      d_dWdx_dW = HALF*(dWnNWdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*(dWnNWdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if (cell_orientation == SOUTH || cell_orientation == SOUTH_WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*(norm[0].y+norm[3].y))/AREA;      
    } /* endif */
    break;
  }
}

/*!
 * Change a block's BC's.
 */
void Grid2D_Quad_Block_HO::set_BCs(const int& FACE, const int& BC){
  switch(FACE){
  case NORTH:
    for ( int i = ICl - Nghost; i <=  ICu + Nghost; i++) {
      BCtypeN[i] = BC;
    } /* endfor */
    break;
  case SOUTH:
    for ( int i = ICl - Nghost; i <=  ICu + Nghost; i++) {
      BCtypeS[i] = BC;
    } /* endfor */
    break;
  case EAST:
    for ( int j = JCl - Nghost; j <=  JCu + Nghost; j++) {
      BCtypeE[j] = BC;
    } /* endfor */
    break;
  case WEST:
    for ( int j = JCl - Nghost; j <=  JCu + Nghost; j++) {
      BCtypeW[j] = BC;
    } /* endfor */
    break;
  default:
    cerr<<"\n WRONG \"FACE\" in set_BCs. \n";
    throw runtime_error("Grid2D_Quad_Block_HO::set_BCs() ERROR! Wrong face!");
  };
}

/*!
 *
 * Create quadrilateral grid block for a Cartesian      
 * mesh defined on a square box with four corner        
 * coordinates (X,Y): (-1,1), (1,1), (1,-1), (-1,1).    
 *                                                      
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block(const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells){
  
  int i, j;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
  Vector2D x_north, x_south, x_east, x_west;
 
  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */
  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);
  
  /* Allocate (re-allocate) memory for the boundary splines 
     defining this quadrilateral mesh block. */
  
  if (BndNorthSpline.np != 0) BndNorthSpline.deallocate();
  BndNorthSpline.allocate(2);
  
  if (BndSouthSpline.np != 0) BndSouthSpline.deallocate();
  BndSouthSpline.allocate(2);
  
  if (BndEastSpline.np != 0) BndEastSpline.deallocate();
  BndEastSpline.allocate(2);
  
  if (BndWestSpline.np != 0) BndWestSpline.deallocate();
  BndWestSpline.allocate(2);
  
  /* Set the boundary spline types to linear, assign
     spline points, and calculate the spline pathlengths. */
  
  BndNorthSpline.settype(SPLINE2D_LINEAR);
  BndNorthSpline.Xp[0] = Vector2D(-ONE, ONE);
  BndNorthSpline.Xp[1] = Vector2D( ONE, ONE);
  BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndNorthSpline.bc[0] = BC_REFLECTION;
  BndNorthSpline.bc[1] = BC_REFLECTION;
  BndNorthSpline.pathlength();
  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];
  
  BndSouthSpline.settype(SPLINE2D_LINEAR);
  BndSouthSpline.Xp[0] = Vector2D(-ONE,-ONE);
  BndSouthSpline.Xp[1] = Vector2D( ONE,-ONE);
  BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndSouthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndSouthSpline.bc[0] = BC_REFLECTION;
  BndSouthSpline.bc[1] = BC_REFLECTION;
  BndSouthSpline.pathlength();
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];
  
  BndEastSpline.settype(SPLINE2D_LINEAR);
  BndEastSpline.Xp[0] = Vector2D( ONE,-ONE);
  BndEastSpline.Xp[1] = Vector2D( ONE, ONE);
  BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.bc[0] = BC_REFLECTION;
  BndEastSpline.bc[1] = BC_REFLECTION;
  BndEastSpline.pathlength();
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];
  
  BndWestSpline.settype(SPLINE2D_LINEAR);
  BndWestSpline.Xp[0] = Vector2D(-ONE,-ONE);
  BndWestSpline.Xp[1] = Vector2D(-ONE, ONE);
  BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.bc[0] = BC_REFLECTION;
  BndWestSpline.bc[1] = BC_REFLECTION;
  BndWestSpline.pathlength();
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];
  
  /* Compute the interior nodes for the quadrilateral mesh block. */
  
  for ( j = JNl ; j <= JNu ; ++j) {
    S_j = double(j-JNl)/double(JNu-JNl);
    
    smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
      BndEastSpline.sp[0];
    s_east = S_j * smax_east + BndEastSpline.sp[0];
    x_east = Spline(s_east, BndEastSpline);
    
    smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
      BndWestSpline.sp[0];
    s_west = S_j * smax_west + BndWestSpline.sp[0];
    x_west = Spline(s_west, BndWestSpline);
    
    for ( i = INl ; i <= INu ; ++i) {
      S_i = double(i-INl)/double(INu-INl);
      
      smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
	BndNorthSpline.sp[0];
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      
      smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
	BndSouthSpline.sp[0];
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);
      
      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;
	
	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;
        
	w_south = (ONE-s_j/smax_j);
	w_north = (s_j/smax_j);
	w_total = w_south + w_north;
	Node[i][j].X = (w_south*x_south + 
			     w_north*x_north)/w_total;
      } /* endif */
      
    } /* endfor */
  }/* endfor */
  
  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();
  
  /* Set the boundary condition types at the grid block
     boundaries. */
  
  Set_BCs();
  
  /* Compute the exterior nodes for the quadrilateral mesh block. */
  
  Update_Exterior_Nodes();
  
  /* Compute the cells for the quadrilateral mesh block. */
  
  Update_Cells();
  
}

/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given in the input boundary spline file:   
 *                                                      
 * Bnd_Spline_File_Name_ptr.                            
 *                                                      
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block(char *Bnd_Spline_File_Name_ptr,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells) {
  
  ifstream bnd_spline_file;
  int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
  int node_init_procedure;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;
  
  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */

  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);

  /* Open the data file containing the boundary splines. */

  bnd_spline_file.open(Bnd_Spline_File_Name_ptr,ios::in);

  /* For each of the north, south, east, and west boundaries
     of this mesh block, read in the number of spline points, 
     allocate memory for the boundary splines, read in the 
     spline points, and finally calculate the spline 
     pathlengths. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndNorthSpline.np != 0) BndNorthSpline.deallocate();
  BndNorthSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndNorthSpline.settype(spline_type);
  bnd_spline_file >> BndNorthSpline;
  BndNorthSpline.pathlength();
  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndSouthSpline.np != 0) BndSouthSpline.deallocate();
  BndSouthSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndSouthSpline.settype(spline_type);
  bnd_spline_file >> BndSouthSpline;
  BndSouthSpline.pathlength();
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndEastSpline.np != 0) BndEastSpline.deallocate();
  BndEastSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndEastSpline.settype(spline_type);
  bnd_spline_file >> BndEastSpline;
  BndEastSpline.pathlength();
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> npts;
  bnd_spline_file.unsetf(ios::skipws);
  if (BndWestSpline.np != 0) BndWestSpline.deallocate();
  BndWestSpline.allocate(npts);
  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> spline_type;
  bnd_spline_file.unsetf(ios::skipws);
  BndWestSpline.settype(spline_type);
  bnd_spline_file >> BndWestSpline;
  BndWestSpline.pathlength();
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

  /* Read the node initialization procedure for this 
     quadrilateral grid block. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> node_init_procedure;
  bnd_spline_file.unsetf(ios::skipws);

  /* Input the type of node distribution stretching 
     functions to be used for each of the coordinate
     directions. Also read in the related stretching
     parameters. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> StretchI >> BetaI >> TauI;
  bnd_spline_file >> StretchJ >> BetaJ >> TauJ;
  bnd_spline_file.unsetf(ios::skipws);

  /* Input the grid orthogonality specifiers for
     the north, south, east, and west boundaries. */

  bnd_spline_file.setf(ios::skipws);
  bnd_spline_file >> OrthogonalN >> OrthogonalS
		  >> OrthogonalE >> OrthogonalW;
  bnd_spline_file.unsetf(ios::skipws);

  /* Close the data file containing the boundary splines. */

  bnd_spline_file.close();

  /* Compute the interior nodes for the quadrilateral mesh block. */

  smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
    BndNorthSpline.sp[0];
  smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
    BndSouthSpline.sp[0];

  smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
    BndEastSpline.sp[0];
  smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
    BndWestSpline.sp[0];

  dS_i = ONE/(TEN*double(INu-INl));
  dS_j = ONE/(TEN*double(JNu-JNl));

  for ( j = JNl ; j <= JNu ; ++j) {
    for ( i = INl ; i <= INu ; ++i) {
      S_j = double(j-JNl)/double(JNu-JNl);
      S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
      s_east = S_j * smax_east + BndEastSpline.sp[0];
      x_east = Spline(s_east, BndEastSpline);
      s_west = S_j * smax_west + BndWestSpline.sp[0];
      x_west = Spline(s_west, BndWestSpline);

      S_i = double(i-INl)/double(INu-INl);
      S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);

      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;

	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	switch(node_init_procedure) {
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			       w_north*x_north)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	  w_west =  (ONE-s_i/smax_i);
	  w_east =  (s_i/smax_i);
	  w_total = w_east + w_west;
	  Node[i][j].X = (w_west*x_west + 
			       w_east*x_east)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	  kx_max = 5;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 7;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(2);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(2);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(kx_max);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S1_j;
	      w_north =  S1_j;
	      w_total = w_north + w_south;
	      S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(kx_max);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S2_j;
	      w_north =  S2_j;
	      w_total = w_north + w_south;
	      S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);

	  if (ky == 0 || ky == ky_max - 2) {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	  } else {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	  kx_max = 7;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 5;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(ky_max);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S1_i;
	      w_east =  S1_i;
	      w_total = w_east + w_west;
	      S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(ky_max);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S2_i;
	      w_east =  S2_i;
	      w_total = w_east + w_west;
	      S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(2);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(2);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);

	  if (kx == 0 || kx == kx_max - 2) {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	  } else {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  //============================================================
	default:
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			       w_north*x_north)/w_total;
	  break;
	} /* endswitch */

      } /* endif */

    } /* endfor */
  }/* endfor */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Set the boundary condition types at the quadrilateral 
     grid block boundaries. */

  Set_BCs();

  /* Compute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Compute the cells for the quadrilateral mesh block. */
  
  Update_Cells();
}


/*!
 * Create a 2D quadrilateral grid block with the four   
 * boundaries of the block defined by blended splines   
 * which are given as input arguments to the routine.   
 *                                                      
 */
void Grid2D_Quad_Block_HO::Create_Quad_Block(Spline2D_HO &Bnd_Spline_North,
					     Spline2D_HO &Bnd_Spline_South,
					     Spline2D_HO &Bnd_Spline_East,
					     Spline2D_HO &Bnd_Spline_West,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells,
					     const int Node_Init_Procedure,
					     const int Stretch_I,
					     const double &Beta_I, 
					     const double &Tau_I,
					     const int Stretch_J,
					     const double &Beta_J,
					     const double &Tau_J,
					     const int Orthogonal_North,
					     const int Orthogonal_South,
					     const int Orthogonal_East,
					     const int Orthogonal_West) {
  
  int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;

  /* Allocate (re-allocate) memory for the cells and nodes 
     of the quadrilateral mesh block. */

  allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);

  /* For each of the north, south, east, and west boundaries
     of this mesh block, assign the boundary splines specified
     by the subroutine input parameters. */

  BndNorthSpline = Bnd_Spline_North;
  BndSouthSpline = Bnd_Spline_South;
  BndEastSpline  = Bnd_Spline_East;
  BndWestSpline  = Bnd_Spline_West;

  SminN = BndNorthSpline.sp[0];
  SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];
  SminS = BndSouthSpline.sp[0];
  SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];
  SminE = BndEastSpline.sp[0];
  SmaxE = BndEastSpline.sp[BndEastSpline.np-1];
  SminW = BndWestSpline.sp[0];
  SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

  /* Assign values for the type of node distribution 
     stretching functions to be used for each of the 
     coordinate directions. Also assign values to the
     related stretching parameters. */

  StretchI = Stretch_I;
  BetaI = Beta_I;
  TauI = Tau_I;
  StretchJ = Stretch_J;
  BetaJ = Beta_J;
  TauJ = Tau_J;

  /* Assign values to the grid orthogonality specifiers for
     the north, south, east, and west boundaries. */

  OrthogonalN = Orthogonal_North;
  OrthogonalS = Orthogonal_South;
  OrthogonalE = Orthogonal_East;
  OrthogonalW = Orthogonal_West;

  /* Compute the interior nodes for the quadrilateral mesh block. */

  smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
    BndNorthSpline.sp[0];
  smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
    BndSouthSpline.sp[0];

  smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
    BndEastSpline.sp[0];
  smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
    BndWestSpline.sp[0];

  for ( j = JNl ; j <= JNu ; ++j) {
    for ( i = INl ; i <= INu ; ++i) {
      S_j = double(j-JNl)/double(JNu-JNl);
      S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
      s_east = S_j * smax_east + BndEastSpline.sp[0];
      x_east = Spline(s_east, BndEastSpline);
      s_west = S_j * smax_west + BndWestSpline.sp[0];
      x_west = Spline(s_west, BndWestSpline);

      S_i = double(i-INl)/double(INu-INl);
      S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
      s_north = S_i * smax_north + BndNorthSpline.sp[0];
      x_north = Spline(s_north, BndNorthSpline);
      s_south = S_i * smax_south + BndSouthSpline.sp[0];
      x_south = Spline(s_south, BndSouthSpline);

      if (j == JNl) {
	Node[i][j].X = x_south;
      } else if (j == JNu) {
	Node[i][j].X = x_north;
      } else if (i == INl) {
	Node[i][j].X = x_west;
      } else if (i == INu) {
	Node[i][j].X = x_east;
      } else {
	s_i = (ONE-S_j)*s_south + S_j*s_north;
	s_j = (ONE-S_i)*s_west + S_i*s_east;

	smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	switch(Node_Init_Procedure) {
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			  w_north*x_north)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	  w_west =  (ONE-s_i/smax_i);
	  w_east =  (s_i/smax_i);
	  w_total = w_east + w_west;
	  Node[i][j].X = (w_west*x_west + 
			  w_east*x_east)/w_total;
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	  kx_max = 5;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 7;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(2);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(2);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(kx_max);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S1_j;
	      w_north =  S1_j;
	      w_total = w_north + w_south;
	      S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(kx_max);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
	      s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
	      x_north = Spline(s_north, BndNorthSpline);

	      s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
	      x_south = Spline(s_south, BndSouthSpline);

	      w_south =  ONE-S2_j;
	      w_north =  S2_j;
	      w_total = w_north + w_south;
	      S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
	      S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);
      
	  if (ky == 0 || ky == ky_max - 2) {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
			    w_north*x_north)/w_total;
	  } else {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
			    w_east*x_east)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	  kx_max = 7;
	  dS_i = ONE/double(kx_max-1);
	  kx = int(floor(S_i/dS_i));
	  S1_i = max(ZERO, double(kx)*dS_i);
	  S2_i = min(ONE, S1_i + dS_i);

	  ky_max = 5;
	  dS_j = ONE/double(ky_max-1);
	  ky = int(floor(S_j/dS_j));
	  S1_j = max(ZERO, double(ky)*dS_j);
	  S2_j = min(ONE, S1_j + dS_j);

	  if (S1_i - TOLER > ZERO) {
	    S1_NS.allocate(ky_max);
	    S1_NS.settype(SPLINE2D_LINEAR);

	    s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	    S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S1_i * smax_north + BndNorthSpline.sp[0],
	      S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S1_i;
	      w_east =  S1_i;
	      w_total = w_east + w_west;
	      S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S1_NS.pathlength();
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } else {
	    S1_NS = BndWestSpline;
	    smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
	      S1_NS.sp[0];
	  } /* endif */

	  if (S2_i + TOLER < ONE) {
	    S2_NS.allocate(ky_max);
	    S2_NS.settype(SPLINE2D_LINEAR);

	    s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	    S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	    S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_north = S2_i * smax_north + BndNorthSpline.sp[0],
	      S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	    S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	    for ( k = 1 ; k <= ky_max-2 ; ++k) {
	      s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
	      x_east = Spline(s_east, BndEastSpline);

	      s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
	      x_west = Spline(s_west, BndWestSpline);

	      w_west =  ONE-S2_i;
	      w_east =  S2_i;
	      w_total = w_east + w_west;
	      S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
	      S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	    } /* endfor */

	    S2_NS.pathlength();
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } else {
	    S2_NS = BndEastSpline;
	    smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
	      S2_NS.sp[0];
	  } /* endif */

	  if (S1_j - TOLER > ZERO) {
	    S1_EW.allocate(2);
	    S1_EW.settype(SPLINE2D_LINEAR);

	    s_west = S1_j * smax_west + BndWestSpline.sp[0];
	    S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S1_j * smax_east + BndEastSpline.sp[0];
	    S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S1_EW.pathlength();
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } else {
	    S1_EW = BndSouthSpline;
	    smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
	      S1_EW.sp[0];
	  } /* endif */

	  if (S2_j + TOLER < ONE) {
	    S2_EW.allocate(2);
	    S2_EW.settype(SPLINE2D_LINEAR);

	    s_west = S2_j * smax_west + BndWestSpline.sp[0];
	    S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	    S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	    s_east = S2_j * smax_east + BndEastSpline.sp[0];
	    S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	    S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	    S2_EW.pathlength();
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } else {
	    S2_EW = BndNorthSpline;
	    smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
	      S2_EW.sp[0];
	  } /* endif */

	  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	  x_north = Spline(s_north, S2_EW);
	  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	  x_south = Spline(s_south, S1_EW);

	  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	  x_east = Spline(s_east, S2_NS);
	  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	  x_west = Spline(s_west, S1_NS);
      
	  if (kx == 0 || kx == kx_max - 2) {
	    w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	    w_east =  (S_i-S1_i)/(S2_i-S1_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
			    w_east*x_east)/w_total;
	  } else {
	    w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	    w_north =  (S_j-S1_j)/(S2_j-S1_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
			    w_north*x_north)/w_total;
	  } /* endif */

	  S1_NS.deallocate(); S2_NS.deallocate();
	  S2_EW.deallocate(); S2_EW.deallocate();
	  break;
	  //============================================================
	default:
	  w_south = (ONE-s_j/smax_j);
	  w_north = (s_j/smax_j);
	  w_total = w_south + w_north;
	  Node[i][j].X = (w_south*x_south + 
			  w_north*x_north)/w_total;
	  break;
	} /* endswitch */

      } /* endif */

    } /* endfor */
  }/* endfor */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Set the boundary condition types at the quadrilateral 
     grid block boundaries. */

  Set_BCs();

  /* Compute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Compute the cells for the quadrilateral mesh block. */

  Update_Cells();

}

/*!
 * Broadcast quadrilateral grid block to all            
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.   
 *                                                      
 */
void Grid2D_Quad_Block_HO::Broadcast_Quad_Block(void) {

#ifdef _MPI_VERSION

    int i, j, ni, nj, ng, mesh_allocated, buffer_size;
    double *buffer;
 
    /* Broadcast the number of cells in each direction. */

    if (CFFC_Primary_MPI_Processor()) {
      ni = NCi;
      nj = NCj;
      ng = Nghost;
      if (Node != NULL) {
         mesh_allocated = 1;
      } else {
	 mesh_allocated = 0;
      } /* endif */ 
    } /* endif */

    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&mesh_allocated, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the cells and nodes of the quadrilateral 
       mesh block as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
       if (NCi != ni || NCj != nj || Nghost != ng) { 
          if (mesh_allocated) allocate(ni-2*ng, nj-2*ng, ng); 
       } /* endif */
    } /* endif */

    /* Broadcast the north, south, east, and west 
       boundary splines. */

    BndNorthSpline.Broadcast_Spline();
    BndSouthSpline.Broadcast_Spline();
    BndEastSpline.Broadcast_Spline();
    BndWestSpline.Broadcast_Spline();

   /* Broadcast min/max pathlengths for boundary splines. */

    MPI::COMM_WORLD.Bcast(&(SminN), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SmaxN), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SminS), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SmaxS), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SminE), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SmaxE), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SminW), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(SmaxW), 1, MPI::DOUBLE, 0);

    /* Broadcast node control parameters. */

    MPI::COMM_WORLD.Bcast(&(StretchI), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(BetaI), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(TauI), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(StretchJ), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(BetaJ), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(TauJ), 1, MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(&(OrthogonalN), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(OrthogonalS), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(OrthogonalE), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(OrthogonalW), 1, MPI::INT, 0);

    /* Broadcast the node locations for grid block. */

    if (mesh_allocated) {
       ni = (INu+Nghost) - (INl-Nghost) + 1;
       nj = (JNu+Nghost) - (JNl-Nghost) + 1;
       buffer = new double[2*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
              for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
 	          buffer[buffer_size] = Node[i][j].X.x;
 	          buffer[buffer_size+1] = Node[i][j].X.y;
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 2*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
              for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
 	          Node[i][j].X.x = buffer[buffer_size];
 	          Node[i][j].X.y = buffer[buffer_size+1];
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       delete []buffer;
       buffer = NULL;

    } /* endif */

    /*  On non-primary MPI processors, set the boundary condition types
        compute the exterior nodes, and compute the cells for the 
        quadrilateral mesh block. */

    if (mesh_allocated && !CFFC_Primary_MPI_Processor()) {
      /* Require update of the whole mesh */
      Schedule_Interior_Mesh_Update();
      Schedule_Ghost_Cells_Update();

      Set_BCs();
      Update_Exterior_Nodes();
      Update_Cells();
    } /* endif */
#endif

}

#ifdef _MPI_VERSION
/*!
 * Broadcast quadrilateral grid block to all processors 
 * associated with the specified communicator from the  
 * specified processor using the MPI broadcast routine. 
 */
void Grid2D_Quad_Block_HO::Broadcast_Quad_Block(MPI::Intracomm &Communicator, 
						const int Source_CPU) {
  
  int Source_Rank = 0;
  int i, j, ni, nj, ng, mesh_allocated, buffer_size;
  double *buffer;
  
  /* Broadcast the number of cells in each direction. */
  
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = NCi;
    nj = NCj;
    ng = Nghost;
    if (Node != NULL) {
      mesh_allocated = 1;
    } else {
      mesh_allocated = 0;
    } /* endif */ 
  } /* endif */

  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&mesh_allocated, 1, MPI::INT, Source_Rank);
  
  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the cells and nodes of the quadrilateral 
     mesh block as necessary. */
  
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (NCi != ni || NCj != nj || Nghost != ng) { 
      if (mesh_allocated) allocate(ni-2*ng, nj-2*ng, ng); 
    } /* endif */
  } /* endif */
  
    /* Broadcast the north, south, east, and west 
       boundary splines. */
  
  Broadcast_Spline(BndNorthSpline, Communicator, Source_CPU);
  Broadcast_Spline(BndSouthSpline, Communicator, Source_CPU);
  Broadcast_Spline(BndEastSpline, Communicator, Source_CPU);
  Broadcast_Spline(BndWestSpline, Communicator, Source_CPU);
  
  /* Broadcast min/max pathlengths for boundary splines. */
  
  Communicator.Bcast(&(SminN), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxN), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminS), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxS), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminE), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxE), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SminW), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(SmaxW), 1, MPI::DOUBLE, Source_Rank);
  
  /* Broadcast node control parameters. */

  Communicator.Bcast(&(StretchI), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(BetaI), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(TauI), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(StretchJ), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(BetaJ), 1, MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(TauJ), 1, MPI::DOUBLE, Source_Rank);

  Communicator.Bcast(&(OrthogonalN), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalS), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalE), 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&(OrthogonalW), 1, MPI::INT, Source_Rank);

  /* Broadcast the node locations for grid block. */

  if (mesh_allocated) {
    ni = (INu+Nghost) - (INl-Nghost) + 1;
    nj = (JNu+Nghost) - (JNl-Nghost) + 1;
    buffer = new double[2*ni*nj];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      buffer_size = 0;
      for (j  = JNl-Nghost ; j <= JNu+Nghost ; ++j ) {
	for ( i = INl-Nghost ; i <= INu+Nghost ; ++i ) {
	  buffer[buffer_size] = Node[i][j].X.x;
	  buffer[buffer_size+1] = Node[i][j].X.y;
	  buffer_size = buffer_size + 2;
	} /* endfor */
      } /* endfor */
    } /* endif */
    
    buffer_size = 2*ni*nj;
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);
    
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      buffer_size = 0;
      for (j  = JNl-Nghost; j <= JNu+Nghost; ++j ) {
	for ( i = INl-Nghost; i <= INu+Nghost; ++i ) {
	  Node[i][j].X.x = buffer[buffer_size];
	  Node[i][j].X.y = buffer[buffer_size+1];
	  buffer_size = buffer_size + 2;
	} /* endfor */
      } /* endfor */
    } /* endif */
    
    delete []buffer; 
    buffer = NULL;
    
  } /* endif */
  
    /*  On non-source MPI processors, set the boundary condition types
        compute the exterior nodes, and compute the cells for the 
        quadrilateral mesh block. */
  
  if (mesh_allocated && 
      !(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    
    /* Require update of the whole mesh */
    Schedule_Interior_Mesh_Update();
    Schedule_Ghost_Cells_Update();
    
    Set_BCs();
    Update_Exterior_Nodes();
    Update_Cells();
  } /* endif */

}
#endif

/*!
 * Smooths quadrilateral grid block using an elliptic   
 * grid generation method (transformed Poisson's        
 * equation method) with a successive over-relaxation   
 * (SOR) Gauss-Seidel solution procedure.  The routine  
 * follows the technique outlined by Sorenson (1980)    
 * and enforces fixed grid spacing and orthogonality at  
 * the bottom (south) and top (north) boundaries of the 
 * quadrilateral block as required.  Modifications have 
 * been introduced to permit the enforcement of grid    
 * orthogonality at the left (west) and right (east)    
 * boundaries as well.                                  
 *                                                      
 */
void Grid2D_Quad_Block_HO::Smooth_Quad_Block(const int Number_of_Iterations) {

  // Run smoother only is Smooth_Quad_Block_Flag is ON
  if (Smooth_Quad_Block_Flag){
  
    /*!
     * \verbatim
     *  Local Variable description:
     * 
     *  xij, yij                  Two-dimensional arrays containing
     *                            the coordinates of the nodes for
     *                            the grid block in physical (x,y) space.
     *
     *  deltas_b, deltas_t,       One-dimensional arrays containing the
     *  deltas_l, deltas_r        node spacing in the normal direction
     *                            between the first two nodes at the upper
     *                            and lower boundaries of the grid block.
     *
     *  xz_b, yz_b, xe_b, ye_b,   One-dimensional arrays containing the
     *  xzz_b, yzz_b, xee_b,      partial derivatives of x and y
     *  yee_b, xze_b, yze_b       coordinates of the physical space
     *                            with respect to the zeta and eta
     *  xz_t, yz_t, xe_t, ye_t,   coordinates of the computational space,
     *  xzz_t, yzz_t, xee_t,      all evaluated at the bottom (south),
     *  yee_t, xze_t, yze_t       top (north), left (west), and right 
     *                            (east) boundaries of the grid block.
     *
     *  xz_l, yz_l, xe_l, ye_l,
     *  xzz_l, yzz_l, xee_l,   
     *  yee_l, xze_l, yze_l
     *
     *  xz_r, yz_r, xe_r, ye_r,
     *  xzz_r, yzz_r, xee_r,   
     *  yee_r, xze_r, yze_r
     *
     *  jacob_b, jacob_t,         One-dimensional arrays containing the
     *  jacob_l, jacob_r          Jacobian of the coordinate transformation
     *                            evaluated at the bottom (south), top 
     *                            (north), left (west), and right (east)
     *                            boundaries of the grid block.
     *
     *  pb, pt, pl, pr,           One-dimensional arrays containing the
     *  qb, qt, ql, qr            values of the source terms of the
     *                            transformed Poisson's equations at the
     *                            bottom (south), top (north), left (west), 
     *                            and right (east) boundaries of the grid 
     *                            block.
     *
     *  norm_b, norm_t,           Integer parameters indicating whether 
     *  norm_l, norm_r            orthogonality at the bottom (south), top
     *                            (north), left (west), and right (east) 
     *                            boundaries is to be enforced.  If value 
     *                            is zero then orthogonality is not enforced 
     *                            and any other value causes the orthogonality 
     *                            condition to be applied.
     *
     *  wsor, wsur                SOR and successive under-relaxation (SUR)
     *                            iteration parameters.  wsor is used to
     *                            over-relax the solution for x and y and
     *                            wsur is used to under-relax the solution
     *                            for p, q, r, and s.  Typically, 
     *                            wsor=1.5-1.8 and wsur=0.30-0.70.
     *
     *  a, b,                     Exponential decay parameters defining the
     *  c, d,                     effects of the source terms away from the
     *  e                         bottom (south), top (north), left (west), 
     *                            and right (east) boundaries, as well as the 
     *                            corners.  If a (or b, c, d, e) is small 
     *                            (i.e., a = 0.2) source terms have effects far 
     *                            from boundary.  If a (or b, c, d, e) is large
     *                            (i.e., a = 0.7) source terms have less 
     *                            effects far from boundary.
     *
     *  n_zeta, n_eta             Number of nodes in the zeta and eta
     *                            directions of the computational domain for 
     *                            grid block.
     *
     * alpha, beta, gamma        Parameters used in the transformed Poisson's
     *                            equations.
     * r1, r2, r3, r4            Parameters used to evaluate the source terms
     *                            of the transformed Poisson's equations at
     *                            the bottom (south), top (north), left (west), 
     *                            and right (east) boundaries of the grid block.
     *
     * dpb, dpt, dpl, dpr,       Changes in pb, pt, pl, pr, qb, qt, ql, and qr
     * dqb, dqt, dql, dqr        before the under-relaxation and limiting
     *                           algorithm is applied.
     *
     * dpmax, dqmax              Maximum permitted changes in p and q at the
     *                           bottom, top, left, and right boundaries.
     *
     * xzij, yzij, xeij, yeij,   Local values of the transformation metrics
     * jacobij                   and Jacobian.
     *
     * pij, qij                  Local values of the source terms.
     *
     * dxij, dyij                Changes in xij and yij for on step in the
     *                           iteration technique before the over-
     *                           relaxation algorithm is applied.
     *
     * fb, ft, fl, fr,           Exponentials use to evaluate the local
     * fbl, fbr, ftl, ftr        values of the source terms. 
     *
     * \endverbatim
     */

    double **xij, **yij,
      *deltas_b, *deltas_t, *deltas_l, *deltas_r,
      *xz_b, *yz_b, *xe_b, *ye_b, 
      *xzz_b, *yzz_b, *xee_b, *yee_b, *xze_b, *yze_b,
      *xz_t, *yz_t, *xe_t, *ye_t,
      *xzz_t, *yzz_t, *xee_t, *yee_t, *xze_t, *yze_t,
      *xz_l, *yz_l, *xe_l, *ye_l,
      *xzz_l, *yzz_l, *xee_l, *yee_l, *xze_l, *yze_l,
      *xz_r, *yz_r, *xe_r, *ye_r,
      *xzz_r, *yzz_r, *xee_r, *yee_r, *xze_r, *yze_r,
      *pb, *pt, *pl, *pr,
      *qb, *qt, *ql, *qr,
      *jacob_b, *jacob_t, *jacob_l, *jacob_r;

    int num_iter, i_poisson, 
      norm_b, norm_t, norm_l, norm_r,
      n_zeta, n_eta;
    double wsor, wsur, a, b, c, d, e, l_damping;

    int i, j, i_step, ii, jj;
    double alpha, beta, gamma;
    double r1, r2, r3, r4;
    double dpb, dpt, dpl, dpr, dqb, dqt, dql, dqr;
    double dpmax, dqmax;
    double xzij, yzij, xeij, yeij, jacobij, pij, qij, dxij, dyij;
    double fb, ft, fl, fr, fbl, fbr, ftl, ftr;

    /* Determine the size of and create the Poisson equation solution 
       variables for the elliptic smoothing algorithm. */

    n_zeta = INu - INl + 1;
    n_eta = JNu - JNl + 1;

    xij = new double*[n_zeta];
    yij = new double*[n_zeta];
    for ( i = 0; i <= n_zeta-1 ; ++i ) {
      xij[i] = new double[n_eta];
      yij[i] = new double[n_eta];
    } /* endfor */

    deltas_b = new double[n_zeta];
    deltas_t = new double[n_zeta];
    deltas_l = new double[n_eta];
    deltas_r = new double[n_eta];

    xz_b = new double[n_zeta]; 
    yz_b = new double[n_zeta];
    xe_b = new double[n_zeta];
    ye_b = new double[n_zeta];
    xzz_b = new double[n_zeta];
    yzz_b = new double[n_zeta];
    xee_b = new double[n_zeta];
    yee_b = new double[n_zeta];
    xze_b = new double[n_zeta];
    yze_b = new double[n_zeta];

    xz_t = new double[n_zeta];
    yz_t = new double[n_zeta];
    xe_t = new double[n_zeta]; 
    ye_t = new double[n_zeta];
    xzz_t = new double[n_zeta];
    yzz_t = new double[n_zeta];
    xee_t = new double[n_zeta];
    yee_t = new double[n_zeta];
    xze_t = new double[n_zeta];
    yze_t = new double[n_zeta];

    xz_l = new double[n_eta];
    yz_l = new double[n_eta];
    xe_l = new double[n_eta];
    ye_l = new double[n_eta];
    xzz_l = new double[n_eta];
    yzz_l = new double[n_eta];
    xee_l = new double[n_eta];
    yee_l = new double[n_eta];
    xze_l = new double[n_eta];
    yze_l = new double[n_eta];

    xz_r = new double[n_eta]; 
    yz_r = new double[n_eta];
    xe_r = new double[n_eta];
    ye_r = new double[n_eta];
    xzz_r = new double[n_eta];
    yzz_r = new double[n_eta];
    xee_r = new double[n_eta];
    yee_r = new double[n_eta];
    xze_r = new double[n_eta];
    yze_r = new double[n_eta];

    pb = new double[n_zeta];
    pt = new double[n_zeta];
    pl = new double[n_eta];
    pr = new double[n_eta];

    qb = new double[n_zeta]; 
    qt = new double[n_zeta]; 
    ql = new double[n_eta];
    qr = new double[n_eta];

    jacob_b = new double[n_zeta];
    jacob_t = new double[n_zeta];
    jacob_l = new double[n_eta];
    jacob_r = new double[n_eta];

    /* Assign initial values to the Poisson equation solution
       variables using the current node locations for the nodes
       in the quadrilateral grid block. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	ii = i - INl;
	jj = j - JNl;
	xij[ii][jj] = Node[i][j].X.x;
	yij[ii][jj] = Node[i][j].X.y;
      } /* endfor */
    } /* endfor */

    /* Set the Poisson equation solution parameters. */

    norm_b = OrthogonalS;
    norm_t = OrthogonalN;
    norm_l = OrthogonalW;
    norm_r = OrthogonalE;

    wsor = 1.25;
    wsur = 0.10;
    l_damping = 1.00;

    if (norm_b == 0 && norm_t == 0 &&
	norm_l == 0 && norm_r == 0) {
      i_poisson = 0;
      a = ONE;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
	       norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
	       norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } /* endif */

    /* Initialize the partial derivatives and source terms at the 
       bottom (south), top (north), left (west), and right (east) 
       boundaries.  The procedure for doing this depends on whether 
       or not orthogonality is enforced at the bottom (south), 
       top (north), left (west), and right (east) boundaries. */

    /* BOTTOM AND TOP */

    xz_b[0] = HALF * (xij[1][0] - xij[0][0]);
    xz_t[0] = HALF * (xij[1][n_eta-1] - xij[0][n_eta-1]);

    yz_b[0] = HALF * (yij[1][0] - yij[0][0]);
    yz_t[0] = HALF * (yij[1][n_eta-1] - yij[0][n_eta-1]);

    xe_b[0] = ZERO;
    xe_t[0] = ZERO;

    ye_b[0] = ZERO;
    ye_t[0] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xz_b[i] = HALF * (xij[i+1][0] - xij[i-1][0]);
      xz_t[i] = HALF * (xij[i+1][n_eta-1] - xij[i-1][n_eta-1]);
 
      yz_b[i] = HALF * (yij[i+1][0] - yij[i-1][0]);
      yz_t[i] = HALF * (yij[i+1][n_eta-1] - yij[i-1][n_eta-1]);

      xzz_b[i] = xij[i+1][0] - TWO * xij[i][0] + xij[i-1][0];
      xzz_t[i] = xij[i+1][n_eta-1] - TWO * xij[i][n_eta-1] +
	xij[i-1][n_eta-1];

      yzz_b[i] = yij[i+1][0] - TWO * yij[i][0] + yij[i-1][0];
      yzz_t[i] = yij[i+1][n_eta-1] - TWO * yij[i][n_eta-1] +
	yij[i-1][n_eta-1];
 
      if (norm_b == 0) {
	xe_b[i] = xij[i][1] - xij[i][0];
	ye_b[i] = yij[i][1] - yij[i][0];
      } else {
	deltas_b[i] = hypot( (xij[i][1]-xij[i][0]),
			     (yij[i][1]-yij[i][0]) );
	xe_b[i] = - deltas_b[i] * yz_b[i] /
	  hypot( xz_b[i], yz_b[i] );
	ye_b[i] = deltas_b[i] * xz_b[i] /
	  hypot( xz_b[i], yz_b[i] );
      } /* endif */

      if (norm_t == 0) {
	xe_t[i] = xij[i][n_eta-1] - xij[i][n_eta-2];
	ye_t[i] = yij[i][n_eta-1] - yij[i][n_eta-2];
      } else {
	deltas_t[i] = hypot( (xij[i][n_eta-1]-xij[i][n_eta-2]),
			     (yij[i][n_eta-1]-yij[i][n_eta-2]) );
	xe_t[i] = - deltas_t[i] * yz_t[i] /
	  hypot( xz_t[i], yz_t[i] );
	ye_t[i] = deltas_t[i] * xz_t[i] /
	  hypot( xz_t[i], yz_t[i] );
      } /* endif */      

      jacob_b[i] = xz_b[i] * ye_b[i] - xe_b[i] * yz_b[i];
      jacob_t[i] = xz_t[i] * ye_t[i] - xe_t[i] * yz_t[i];

      pb[i] = ZERO;
      pt[i] = ZERO;

      qb[i] = ZERO;
      qt[i] = ZERO;
    } /* endfor */

    xz_b[n_zeta-1] = HALF * (xij[n_zeta-1][0] - xij[n_zeta-2][0]);
    xz_t[n_zeta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
			     xij[n_zeta-2][n_eta-1]);

    yz_b[n_zeta-1] = HALF * (yij[n_zeta-1][0] - yij[n_zeta-2][0]);
    yz_t[n_zeta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
			     yij[n_zeta-2][n_eta-1]);

    xe_b[n_zeta-1] = ZERO;
    xe_t[n_zeta-1] = ZERO;

    ye_b[n_zeta-1] = ZERO;
    ye_t[n_zeta-1] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xze_b[i] = HALF * (xe_b[i+1] - xe_b[i-1]);
      xze_t[i] = HALF * (xe_t[i+1] - xe_t[i-1]);
 
      yze_b[i] = HALF * (ye_b[i+1] - ye_b[i-1]);
      yze_t[i] = HALF * (ye_t[i+1] - ye_t[i-1]);
    } /* endfor */

    /* LEFT AND RIGHT */

    xe_l[0] = HALF * (xij[0][1] - xij[0][0]);
    xe_r[0] = HALF * (xij[n_zeta-1][1] - xij[n_zeta-1][0]);

    ye_l[0] = HALF * (yij[0][1] - yij[0][0]);
    ye_r[0] = HALF * (yij[n_zeta-1][1] - yij[n_zeta-1][0]);

    xz_l[0] = ZERO;
    xz_r[0] = ZERO;

    yz_l[0] = ZERO;
    yz_r[0] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
      xe_l[j] = HALF * (xij[0][j+1] - xij[0][j-1]);
      xe_r[j] = HALF * (xij[n_zeta-1][j+1] - xij[n_zeta-1][j-1]);

      ye_l[j] = HALF * (yij[0][j+1] - yij[0][j-1]);
      ye_r[j] = HALF * (yij[n_zeta-1][j+1] - yij[n_zeta-1][j-1]);

      xee_l[j] = xij[0][j+1] - TWO * xij[0][j] + xij[0][j-1];
      xee_r[j] = xij[n_zeta-1][j+1] - TWO * xij[n_zeta-1][j] +
	xij[n_zeta-1][j-1];

      yee_l[j] = yij[0][j+1] - TWO * yij[0][j] + yij[0][j-1];
      yee_r[j] = yij[n_zeta-1][j+1] - TWO * yij[n_zeta-1][j] +
	yij[n_zeta-1][j-1];

      if (norm_l == 0) {
	xz_l[j] = xij[1][j] - xij[0][j];
	yz_l[j] = yij[1][j] - yij[0][j];
      } else {
	deltas_l[j] = hypot( (xij[1][j]-xij[0][j]),
			     (yij[1][j]-yij[0][j]) );
	xz_l[j] = deltas_l[j] * ye_l[j] /
	  hypot( xe_l[j], ye_l[j] );
	yz_l[j] = - deltas_l[j] * xe_l[j] /
	  hypot( xe_l[j], ye_l[j] );
      } /* endif */

      if (norm_r == 0) {
	xz_r[j] = xij[n_zeta-1][j] - xij[n_zeta-2][j];
	yz_r[j] = yij[n_zeta-1][j] - yij[n_zeta-2][j];
      } else {
	deltas_r[j] = hypot( (xij[n_zeta-1][j]-xij[n_zeta-2][j]),
			     (yij[n_zeta-1][j]-yij[n_zeta-2][j]) );
	xz_r[j] = deltas_r[j] * ye_r[j] /
	  hypot( xe_r[j], ye_r[j] );
	yz_r[j] = -deltas_r[j] * xe_r[j] /
	  hypot( xe_r[j], ye_r[j] );
      } /* endif */

      jacob_l[j] = xz_l[j] * ye_l[j] - xe_l[j] * yz_l[j];
      jacob_r[j] = xz_r[j] * ye_r[j] - xe_r[j] * yz_r[j];

      pl[j] = ZERO;
      pr[j] = ZERO;

      ql[j] = ZERO;
      qr[j] = ZERO;
    } /* endfor */

    xe_l[n_eta-1] = HALF * (xij[0][n_eta-1] - xij[0][n_eta-2]);
    xe_r[n_eta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
			    xij[n_zeta-1][n_eta-2]);

    ye_l[n_eta-1] = HALF * (yij[0][n_eta-1] - yij[0][n_eta-2]);
    ye_r[n_eta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
			    yij[n_zeta-1][n_eta-2]);

    xz_l[n_eta-1] = ZERO;
    xz_r[n_eta-1] = ZERO;

    yz_l[n_eta-1] = ZERO;
    yz_r[n_eta-1] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
      xze_l[j] = HALF * (xz_l[j+1] - xz_l[j-1]);
      xze_r[j] = HALF * (xz_r[j+1] - xz_r[j-1]);

      yze_l[j] = HALF * (yz_l[j+1] - yz_l[j-1]);
      yze_r[j] = HALF * (yz_r[j+1] - yz_r[j-1]);
    } /* endfor */

    /* Begin a new iteration cycle and use the current values of xij
       and yij to determine the derivatives xee and yee at the bottom
       (south) and top (north) boundaries as well as the derivatives
       xzz and yzz at the left (west) and right (east) boundaries, 
       repectively.  These values may then be used to determine pb, 
       pt, pl, pr, qb, qt, ql, and qr.  */

    num_iter = 0;

  next_iteration: ;

    num_iter = num_iter + 1;

    if (i_poisson == 0) goto no_source_term_boundary_evaluation;

    /* BOTTOM AND TOP */

    for ( i = 1; i < n_zeta - 1; ++i ) {
      xee_b[i] = HALF * (-SEVEN*xij[i][0] + EIGHT*xij[i][1] -
			 xij[i][2]) - THREE * xe_b[i];
      xee_t[i] = HALF * (-SEVEN*xij[i][n_eta-1] + 
			 EIGHT*xij[i][n_eta-2] -
			 xij[i][n_eta-3]) + THREE * xe_t[i];

      yee_b[i] = HALF * (-SEVEN*yij[i][0] + EIGHT*yij[i][1] -
			 yij[i][2]) - THREE * ye_b[i];
      yee_t[i] = HALF * (-SEVEN*yij[i][n_eta-1] + 
			 EIGHT*yij[i][n_eta-2] -
			 yij[i][n_eta-3]) + THREE * ye_t[i];

      alpha = xe_b[i] * xe_b[i] + ye_b[i] * ye_b[i];
      beta  = xz_b[i] * xe_b[i] + yz_b[i] * ye_b[i];
      gamma = xz_b[i] * xz_b[i] + yz_b[i] * yz_b[i];
      r1 = - ONE * (alpha * xzz_b[i] - TWO * beta * xze_b[i] +
		    gamma * xee_b[i]) / (jacob_b[i] * jacob_b[i]);
      r2 = - ONE * (alpha * yzz_b[i] - TWO * beta * yze_b[i] +
		    gamma * yee_b[i]) / (jacob_b[i] * jacob_b[i]);

      dpb = (ye_b[i] * r1 - xe_b[i] * r2) / jacob_b[i] - pb[i];
      if (fabs(pb[i]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pb[i]);
      } /* endif */
      if (wsur * fabs(dpb) <= dpmax) {
	pb[i] = pb[i] + wsur * dpb;
      } else {
	if (dpb >= ZERO) {
	  pb[i] = pb[i] + dpmax;
	} else {
	  pb[i] = pb[i] - dpmax;
	} /* endif */
      } /* endif */

      dqb = (xz_b[i] * r2 - yz_b[i] * r1) / jacob_b[i] - qb[i];
      if (fabs(qb[i]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qb[i]);
      } /* endif */
      if (wsur * fabs(dqb) <= dqmax) {
	qb[i] = qb[i] + wsur * dqb;
      } else {
	if (dqb >= ZERO) {
	  qb[i] = qb[i] + dqmax;
	} else {
	  qb[i] = qb[i] - dqmax;
	} /* endif */
      } /* endif */

      alpha = xe_t[i] * xe_t[i] + ye_t[i] * ye_t[i];
      beta  = xz_t[i] * xe_t[i] + yz_t[i] * ye_t[i];
      gamma = xz_t[i] * xz_t[i] + yz_t[i] * yz_t[i];
      r3 = - ONE * (alpha * xzz_t[i] - TWO * beta * xze_t[i] +
		    gamma * xee_t[i]) / (jacob_t[i] * jacob_t[i]);
      r4 = - ONE * (alpha * yzz_t[i] - TWO * beta * yze_t[i] +
		    gamma * yee_t[i]) / (jacob_t[i] * jacob_t[i]);

      dpt = (ye_t[i] * r3 - xe_t[i] * r4) / jacob_t[i] - pt[i];
      if (fabs(pt[i]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pt[i]);
      } /* endif */
      if (wsur * fabs(dpt) <= dpmax) {
	pt[i] = pt[i] + wsur * dpt;
      } else {
	if (dpt >= ZERO) {
	  pt[i] = pt[i] + dpmax;
	} else {
	  pt[i] = pt[i] - dpmax;
	} /* endif */
      } /* endif */

      dqt = (xz_t[i] * r4 - yz_t[i] * r3) / jacob_t[i] - qt[i];
      if (fabs(qt[i]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qt[i]);
      } /* endif */
      if (wsur * fabs(dqt) <= dqmax) {
	qt[i] = qt[i] + wsur * dqt;
      } else {
	if (dqt >= ZERO) {
	  qt[i] = qt[i] + dqmax;
	} else {
	  qt[i] = qt[i] - dqmax;
	} /* endif */
      } /* endif */
    } /* endfor */

    /* LEFT AND RIGHT */

    for ( j = 1; j < n_eta - 1; ++j ) {
      xzz_l[j] = HALF * (-SEVEN*xij[0][j] + EIGHT*xij[1][j] -
			 xij[2][j]) - THREE * xz_l[j];
      xzz_r[j] = HALF * (-SEVEN*xij[n_zeta-1][j] + 
			 EIGHT*xij[n_zeta-2][j] -
			 xij[n_zeta-3][j]) + THREE * xz_r[j];

      yzz_l[j] = HALF * (-SEVEN*yij[0][j] + EIGHT*yij[1][j] -
			 yij[2][j]) - THREE * yz_l[j];
      yzz_r[j] = HALF * (-SEVEN*yij[n_zeta-1][j] + 
			 EIGHT*yij[n_zeta-2][j] -
			 yij[n_zeta-3][j]) + THREE * yz_r[j];

      alpha = xe_l[j] * xe_l[j] + ye_l[j] * ye_l[j];
      beta  = xz_l[j] * xe_l[j] + yz_l[j] * ye_l[j];
      gamma = xz_l[j] * xz_l[j] + yz_l[j] * yz_l[j];
      r1 = - ONE * (alpha * xzz_l[j] - TWO * beta * xze_l[j] +
		    gamma * xee_l[j]) / (jacob_l[j] * jacob_l[j]);
      r2 = - ONE * (alpha * yzz_l[j] - TWO * beta * yze_l[j] +
		    gamma * yee_l[j]) / (jacob_l[j] * jacob_l[j]);

      dpl = (ye_l[j] * r1 - xe_l[j] * r2) / jacob_l[j] - pl[j];
      if (fabs(pl[j]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pl[j]);
      } /* endif */
      if (wsur * fabs(dpl) <= dpmax) {
	pl[j] = pl[j] + wsur * dpl;
      } else {
	if (dpl >= ZERO) {
	  pl[j] = pl[j] + dpmax;
	} else {
	  pl[j] = pl[j] - dpmax;
	} /* endif */
      } /* endif */

      dql = (xz_l[j] * r2 - yz_l[j] * r1) / jacob_l[j] - ql[j];
      if (fabs(ql[j]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(ql[j]);
      } /* endif */
      if (wsur * fabs(dql) <= dqmax) {
	ql[j] = ql[j] + wsur * dql;
      } else {
	if (dql >= ZERO) {
	  ql[j] = ql[j] + dqmax;
	} else {
	  ql[j] = ql[j] - dqmax;
	} /* endif */
      } /* endif */

      alpha = xe_r[j] * xe_r[j] + ye_r[j] * ye_r[j];
      beta  = xz_r[j] * xe_r[j] + yz_r[j] * ye_r[j];
      gamma = xz_r[j] * xz_r[j] + yz_r[j] * yz_r[j];
      r3 = - ONE * (alpha * xzz_r[j] - TWO * beta * xze_r[j] +
		    gamma * xee_r[j]) / (jacob_r[j] * jacob_r[j]);
      r4 = - ONE * (alpha * yzz_r[j] - TWO * beta * yze_r[j] +
		    gamma * yee_r[j]) / (jacob_r[j] * jacob_r[j]);

      dpr = (ye_r[j] * r3 - xe_r[j] * r4) / jacob_r[j] - pr[j];
      if (fabs(pr[j]) <= ONE ) {
	dpmax = HALF;
      } else {
	dpmax = HALF * fabs(pr[j]);
      } /* endif */
      if (wsur * fabs(dpr) <= dpmax) {
	pr[j] = pr[j] + wsur * dpr;
      } else {
	if (dpr >= ZERO) {
	  pr[j] = pr[j] + dpmax;
	} else {
	  pr[j] = pr[j] - dpmax;
	} /* endif */
      } /* endif */

      dqr = (xz_r[j] * r4 - yz_r[j] * r3) / jacob_r[j] - qr[j];
      if (fabs(qr[j]) <= ONE ) {
	dqmax = HALF;
      } else {
	dqmax = HALF * fabs(qr[j]);
      } /* endif */
      if (wsur * fabs(dqr) <= dqmax) {
	qr[j] = qr[j] + wsur * dqr;
      } else {
	if (dqr >= ZERO) {
	  qr[j] = qr[j] + dqmax;
	} else {
	  qr[j] = qr[j] - dqmax;
	} /* endif */
      } /* endif */
    } /* endfor */

    /* Perform one step in the SOR Gauss-Seidel centred 
       finite-difference solution technique. */

  no_source_term_boundary_evaluation: ;

    for ( j = 1; j < n_eta - 1; ++j ) {
      if (i_poisson != 0) {
	fb = exp ( - a * double(j) );
	ft = exp ( b * (double(j) - double(n_eta) + ONE) );
      } /* endfor */

      for ( i = 1; i < n_zeta - 1; ++i ) { 
	if (i_poisson != 0) {
	  fl = exp ( - c * double(i) );
	  fr = exp ( d * (double(i) - double(n_zeta) + ONE) );
 
	  fbl = ONE -
	    exp ( -e * sqrt( double(i)*double(i) +
			     double(j)*double(j) ) );
	  fbr = ONE -
	    exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
			     (double(i) - double(n_zeta) + ONE) +
			     double(j)*double(j) ) );
	  ftl = ONE -
	    exp ( -e * sqrt( double(i)*double(i) +
			     (double(j) - double(n_eta) + ONE) *
			     (double(j) - double(n_eta) + ONE) ) );
	  ftr = ONE -
	    exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
			     (double(i) - double(n_zeta) + ONE) +
			     (double(j) - double(n_eta) + ONE) *
			     (double(j) - double(n_eta) + ONE) ) );
	  pij = (pb[i] * fb + pt[i] * ft + pl[j] * fl + pr[j] * fr) *
	    (fbl * fbr * ftl * ftr);
	  qij = (qb[i] * fb + qt[i] * ft + ql[j] * fl + qr[j] * fr) *
	    (fbl * fbr * ftl * ftr);
	} else {
	  pij = ZERO;
	  qij = ZERO;
	} /* endif */

	xzij = HALF * (xij[i+1][j] - xij[i-1][j]);
	yzij = HALF * (yij[i+1][j] - yij[i-1][j]);
	xeij = HALF * (xij[i][j+1] - xij[i][j-1]);
	yeij = HALF * (yij[i][j+1] - yij[i][j-1]);

	jacobij = xzij * yeij - xeij * yzij;

	alpha = xeij * xeij + yeij * yeij;
	beta  = xzij * xeij + yzij * yeij;
	gamma = xzij * xzij + yzij * yzij;

	dxij = HALF * (jacobij * jacobij * (xzij * pij +
					    xeij * qij) + alpha * (xij[i+1][j] + xij[i-1][j]) +
		       gamma * (xij[i][j+1] + xij[i][j-1]) - HALF * beta *
		       (xij[i+1][j+1] - xij[i+1][j-1] - xij[i-1][j+1] +
			xij[i-1][j-1])) / (alpha + gamma) - xij[i][j];
	xij[i][j] = xij[i][j] + wsor * dxij;
          
	dyij = HALF * (jacobij * jacobij * (yzij * pij +
					    yeij * qij) + alpha * (yij[i+1][j] + yij[i-1][j]) +
		       gamma * (yij[i][j+1] + yij[i][j-1]) - HALF * beta *
		       (yij[i+1][j+1] - yij[i+1][j-1] - yij[i-1][j+1] +
			yij[i-1][j-1])) / (alpha + gamma) - yij[i][j];
	yij[i][j] = yij[i][j] + wsor * dyij;
      } /* endfor */
    } /* endfor */

    /* Check to see if the grid block smoothing is complete.
       If not go to next_iteration. */

    if (num_iter < Number_of_Iterations) goto next_iteration;

    /* Save the newly computed interior node locations for 
       the quadrilateral grid block. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	ii = i - INl;
	jj = j - JNl;
	Node[i][j].X.x = xij[ii][jj];
	Node[i][j].X.y = yij[ii][jj];
      } /* endfor */
    }/* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Re-compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Re-compute the cell values for the quadrilateral mesh block. */

    Update_Cells();

    /* Delete (deallocate) the Poisson equation solution variables. */

    for ( i = 0; i <= n_zeta-1 ; ++i ) {
      delete []xij[i]; xij[i] = NULL;
      delete []yij[i]; yij[i] = NULL;
    } /* endfor */
    delete []xij; xij = NULL;
    delete []yij; yij = NULL;

    delete []deltas_b; deltas_b = NULL;
    delete []deltas_t; deltas_t = NULL;
    delete []deltas_l; deltas_l = NULL; 
    delete []deltas_r; deltas_r = NULL;

    delete []xz_b; xz_b = NULL; 
    delete []yz_b; yz_b = NULL;
    delete []xe_b; xe_b = NULL;
    delete []ye_b; ye_b = NULL;
    delete []xzz_b; xzz_b = NULL;
    delete []yzz_b; yzz_b = NULL;
    delete []xee_b; xee_b = NULL;
    delete []yee_b; yee_b = NULL;
    delete []xze_b; xze_b = NULL;
    delete []yze_b; yze_b = NULL;

    delete []xz_t; xz_t = NULL;
    delete []yz_t; yz_t = NULL;
    delete []xe_t; xe_t = NULL; 
    delete []ye_t; ye_t = NULL;
    delete []xzz_t; xzz_t = NULL;
    delete []yzz_t; yzz_t = NULL;
    delete []xee_t; xee_t = NULL;
    delete []yee_t; yee_t = NULL;
    delete []xze_t; xze_t = NULL;
    delete []yze_t; yze_t = NULL;

    delete []xz_l; xz_l = NULL;
    delete []yz_l; yz_l = NULL;
    delete []xe_l; xe_l = NULL;
    delete []ye_l; ye_l = NULL;
    delete []xzz_l; xzz_l = NULL;
    delete []yzz_l; yzz_l = NULL;
    delete []xee_l; xee_l = NULL;
    delete []yee_l; yee_l = NULL;
    delete []xze_l; xze_l = NULL;
    delete []yze_l; yze_l = NULL;

    delete []xz_r; xz_r = NULL;
    delete []yz_r; yz_r = NULL;
    delete []xe_r; xe_r = NULL;
    delete []ye_r; ye_r = NULL;
    delete []xzz_r; xzz_r = NULL;
    delete []yzz_r; yzz_r = NULL;
    delete []xee_r; xee_r = NULL;
    delete []yee_r; yee_r = NULL;
    delete []xze_r; xze_r = NULL;
    delete []yze_r; yze_r = NULL;

    delete []pb; pb = NULL;
    delete []pt; pt = NULL;
    delete []pl; pl = NULL;
    delete []pr; pr = NULL;

    delete []qb; qb = NULL;
    delete []qt; qt = NULL;
    delete []ql; ql = NULL;
    delete []qr; qr = NULL;

    delete []jacob_b; jacob_b = NULL;
    delete []jacob_t; jacob_t = NULL;
    delete []jacob_l; jacob_l = NULL;
    delete []jacob_r; jacob_r = NULL;

  }
}


/**********************************************************************
 * Routine: Smooth_Rocket_Motor                                       *
 **********************************************************************/
void Grid2D_Quad_Block_HO::Smooth_Rocket_Motor(const double &Length_Chamber,
					       const double &Radius_Chamber,
					       const double &Length_Chamber_To_Throat,
					       const double &Length_Nozzle,
					       const double &Radius_Nozzle_Exit,
					       const double &Radius_Nozzle_Throat,
					       const double &Radius_Grain,
					       const int &Nozzle_Type,
					       const double &Stretching_Factor_Idir,
					       const double &Stretching_Factor_Jdir,
					       const int &sector,
					       const int &level,
					       const int &di,const int &dj,
					       const int &ri,const int &rj,
					       const int &NRi,const int &NRj) {
  
  // Exit immediately if the current block is in the chamber.
  if (Node[INl+5][JNl+5].X.x < ZERO) return ;
  
  Smooth_Quad_Block(min(250,4*max(NCi-4,NCi-4)));
  
  return ;
}

/*!
 * Set the boundary condition type for the north,       
 * south, east, and west boundaries of the              
 * quadrilateral mesh block.                            
 *                                                      
 */

void Grid2D_Quad_Block_HO::Set_BCs(void) {

    int i, j, bc_type_left, bc_type_right;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west;

    if (BndNorthSpline.np == 0) {
       for ( i = ICl-Nghost; i <= ICu+Nghost; ++i) {
	   BCtypeN[i] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeN[ICl-1] = BCtype(SminN,
			       BndNorthSpline);
       BCtypeN[ICu+1] = BCtype(SmaxN,
			       BndNorthSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeN[ICl-GCell] = BCtypeN[ICl-GCell+1];
	 BCtypeN[ICu+GCell] = BCtypeN[ICu+GCell-1];
       }
       for ( i = INl ; i <= INu-1 ; ++i) {
	 s_north = getS(Node[i][JNu].X, BndNorthSpline);
	 bc_type_left = BCtype(s_north, BndNorthSpline);
	 s_north = getS(Node[i+1][JNu].X, BndNorthSpline);
	 bc_type_right = BCtype(s_north, BndNorthSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeN[i] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (BndSouthSpline.np == 0) {
       for ( i = ICl-Nghost; i <= ICu+Nghost ; ++i) {
	   BCtypeS[i] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeS[ICl-1] = BCtype(SminS, 
			       BndSouthSpline);
       BCtypeS[ICu+1] = BCtype(SmaxS, 
			       BndSouthSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeS[ICl-GCell] = BCtypeS[ICl-GCell+1];
	 BCtypeS[ICu+GCell] = BCtypeS[ICu+GCell-1];
       }
       for ( i = INl ; i <= INu-1 ; ++i) {
	 s_south = getS(Node[i][JNl].X, BndSouthSpline);
	 bc_type_left = BCtype(s_south, BndSouthSpline);
	 s_south = getS(Node[i+1][JNl].X, BndSouthSpline);
	 bc_type_right = BCtype(s_south, BndSouthSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeS[i] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (BndEastSpline.np == 0) {
       for ( j = JCl-Nghost; j <= JCu+Nghost ; ++j) {
	   BCtypeE[j] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeE[JCl-1] = BCtype(SminE,
			       BndEastSpline);
       BCtypeE[JCu+1] = BCtype(SmaxE, 
			       BndEastSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeE[JCl-GCell] = BCtypeE[JCl-GCell+1];
	 BCtypeE[JCu+GCell] = BCtypeE[JCu+GCell-1];
       }
       for ( j = JNl ; j <= JNu-1 ; ++j) {
	 s_east = getS(Node[INu][j].X, BndEastSpline);
	 bc_type_left = BCtype(s_east, BndEastSpline);
	 s_east = getS(Node[INu][j+1].X, BndEastSpline);
	 bc_type_right = BCtype(s_east, BndEastSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeE[j] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (BndWestSpline.np == 0) {
       for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
	   BCtypeW[j] = BC_NONE;
       } /* endfor */
    } else {
       BCtypeW[JCl-1] = BCtype(SminW, 
			       BndWestSpline);
       BCtypeW[JCu+1] = BCtype(SmaxW, 
			       BndWestSpline);
       for (int GCell=2; GCell<= Nghost; ++GCell){
	 BCtypeW[JCl-GCell] = BCtypeW[JCl-GCell+1];
	 BCtypeW[JCu+GCell] = BCtypeW[JCu+GCell-1];
       }
       for ( j = JNl ; j <= JNu-1 ; ++j) {
	 s_west = getS(Node[INl][j].X, BndWestSpline);
	 bc_type_left = BCtype(s_west, BndWestSpline);
	 s_west = getS(Node[INl][j+1].X, BndWestSpline);
	 bc_type_right = BCtype(s_west, BndWestSpline);

	 if (bc_type_left == bc_type_right) {
	   BCtypeW[j] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

}

/*!
 * Updates the exterior nodes for the quadrilateral     
 * mesh block.                                          
 */
void Grid2D_Quad_Block_HO::Update_Exterior_Nodes(void) {

    int i, j;
    Vector2D norm_dir, X_norm, X_tan;

    // Update West and East boundary nodes.
    for ( j = JNl ; j <= JNu ; ++j) {
      if (BCtypeW[j] == BCtypeW[j-1] &&
	  BCtypeW[j] != BC_REFLECTION &&
	  BCtypeW[j] != BC_PERIODIC &&
	  BCtypeW[j] != BC_NO_SLIP &&
	  BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j] != BC_MOVING_WALL &&
	  BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j] != BC_BURNING_SURFACE &&
	  BCtypeW[j] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INl+GCell][j].X - 
	     Node[INl][j].X);
	}
      } else if (BCtypeW[j] == BCtypeW[j-1] &&
		 (BCtypeW[j] == BC_REFLECTION ||
		  BCtypeW[j] == BC_NO_SLIP ||
		  BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  BCtypeW[j] == BC_MOVING_WALL ||
		  BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		  BCtypeW[j] == BC_BURNING_SURFACE ||
		  BCtypeW[j] == BC_MASS_INJECTION)) {
	if (j > JNl && j < JNu) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    X_norm = ((Node[INl+GCell][j].X - 
		       Node[INl][j].X) * norm_dir) * norm_dir;
	    X_tan = (Node[INl+GCell][j].X - 
		     Node[INl][j].X) - X_norm;
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      X_norm + X_tan;
	  }
	} else if ( j == JNl ) {
	  //norm_dir = - nfaceW(ICl, j);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      (Node[INl+GCell][j].X - 
	       Node[INl][j].X);
	  }
	} else {
	  //norm_dir = - nfaceW(ICl, j-1);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INl-GCell][j].X = Node[INl][j].X -
	      (Node[INl+GCell][j].X - 
	       Node[INl][j].X);
	  }
	} /* endif */
      } else if (BCtypeW[j] == BCtypeW[j-1] &&
		 BCtypeW[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INu][j].X - 
	     Node[INu-GCell][j].X);
	}
      } /* endif */

      if (BCtypeE[j] == BCtypeE[j-1] &&
	  BCtypeE[j] != BC_REFLECTION &&
	  BCtypeE[j] != BC_PERIODIC &&
	  BCtypeE[j] != BC_NO_SLIP &&
	  BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j] != BC_MOVING_WALL &&
	  BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j] != BC_BURNING_SURFACE &&
	  BCtypeE[j] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j] == BCtypeE[j-1] &&
		 (BCtypeE[j] == BC_REFLECTION ||	  
		  BCtypeE[j] == BC_NO_SLIP ||
		  BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		  BCtypeE[j] == BC_MOVING_WALL ||
		  BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		  BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		  BCtypeE[j] == BC_BURNING_SURFACE ||
		  BCtypeE[j] == BC_MASS_INJECTION)) {
	if (j > JNl && j < JNu) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    X_norm = ((Node[INu][j].X - 
		       Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	    X_tan = (Node[INu][j].X - 
		     Node[INu-GCell][j].X) - X_norm;
	    Node[INu+GCell][j].X = Node[INu][j].X +
	      X_norm - X_tan;
	  }
	} else if ( j == JNl ) {
	  //norm_dir = nfaceE(ICu, j);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INu+GCell][j].X = Node[INu][j].X +
	      (Node[INu][j].X - 
	       Node[INu-GCell][j].X);
	  }
	} else {
	  //norm_dir = nfaceE(ICu, j-1);
	  for(int GCell=1; GCell<=Nghost; ++GCell){
	    Node[INu+GCell][j].X = ( Node[INu][j].X +
				     (Node[INu][j].X - 
				      Node[INu-GCell][j].X) );
	  }
	} /* endif */
      } else if (BCtypeE[j] == BCtypeE[j-1] &&
		 BCtypeE[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */
    
    // Update South and North boundary nodes.
    for ( i = INl ; i <= INu ; ++i) {
	  if (BCtypeS[i] == BCtypeS[i-1] &&
	      BCtypeS[i] != BC_REFLECTION &&
	      BCtypeS[i] != BC_PERIODIC &&
	      BCtypeS[i] != BC_NO_SLIP &&
	      BCtypeS[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[i] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[i] != BC_MOVING_WALL &&
	      BCtypeS[i] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[i] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[i] != BC_BURNING_SURFACE &&
	      BCtypeS[i] != BC_MASS_INJECTION &&
	      BCtypeS[i] != BC_RINGLEB_FLOW) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNl-GCell].X = ( Node[i][JNl].X -
				       (Node[i][JNl+GCell].X - 
					Node[i][JNl].X) );
	    }
	  } else if (BCtypeS[i] == BCtypeS[i-1] &&
		     (BCtypeS[i] == BC_REFLECTION ||  
		      BCtypeS[i] == BC_NO_SLIP ||
		      BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
		      BCtypeS[i] == BC_MOVING_WALL ||
		      BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
		      BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
		      BCtypeS[i] == BC_BURNING_SURFACE ||
		      BCtypeS[i] == BC_MASS_INJECTION ||
		      BCtypeS[i] == BC_RINGLEB_FLOW)) {
	    if (i > INl && i < INu) {
	      //  	      norm_dir = - HALF*(nfaceS(i, JCl) + 
	      //                                  nfaceS(i-1, JCl));
	      // 	      for(int GCell=1; GCell<=Nghost; ++GCell){
	      // 		X_norm = ((Node[i][JNl+GCell].X - 
	      // 			   Node[i][JNl].X) * norm_dir) * norm_dir;
	      // 		X_tan = (Node[i][JNl+GCell].X - 
	      // 			 Node[i][JNl].X) - X_norm;
	      // 		Node[i][JNl-GCell].X = Node[i][JNl].X -
	      // 	                                         X_norm + X_tan;
	      // 	      }
	      if (lfaceS(i,JCl) > NANO &&
		  lfaceS(i-1,JCl) > NANO) {
		norm_dir = - HALF*(nfaceS(i, JCl) + 
				   nfaceS(i-1, JCl));
	      } else if (lfaceS(i,JCl) > NANO) {
		norm_dir = - nfaceS(i, JCl);
	      } else if (lfaceS(i,JCl) > NANO) {
		norm_dir = - nfaceS(i-1, JCl);
	      } else {
	      }
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		X_norm = ((Node[i][JNl+GCell].X - 
			   Node[i][JNl].X) * norm_dir) * norm_dir;
		X_tan = (Node[i][JNl+GCell].X - 
			 Node[i][JNl].X) - X_norm;
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 X_norm + X_tan );
	      }
	    } else if ( i == INl ) {
	      //norm_dir = - nfaceS(i, JCl);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 (Node[i][JNl+GCell].X - 
					  Node[i][JNl].X) );
	      }
	    } else {
	      //norm_dir = - nfaceS(i-1, JCl);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNl-GCell].X = ( Node[i][JNl].X -
					 (Node[i][JNl+GCell].X - 
					  Node[i][JNl].X) );
	      }
	    } /* endif */
	  } else if (BCtypeS[i] == BCtypeS[i-1] &&
		     BCtypeS[i] == BC_PERIODIC) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNl-GCell].X = ( Node[i][JNl].X -
				       (Node[i][JNu].X - 
					Node[i][JNu-GCell].X) );
	  }
        } /* endif */

	  if (BCtypeN[i] == BCtypeN[i-1] &&
	      BCtypeN[i] != BC_REFLECTION &&
	      BCtypeN[i] != BC_PERIODIC &&
	      BCtypeN[i] != BC_NO_SLIP &&
	      BCtypeN[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[i] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[i] != BC_MOVING_WALL &&
	      BCtypeN[i] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[i] != BC_MOVING_WALL_HEATFLUX  &&
	      BCtypeN[i] != BC_BURNING_SURFACE &&
	      BCtypeN[i] != BC_MASS_INJECTION) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNu+GCell].X = ( Node[i][JNu].X +
				       (Node[i][JNu].X - 
					Node[i][JNu-GCell].X) );
	    }
	  } else if (BCtypeN[i] == BCtypeN[i-1] &&
		     (BCtypeN[i] == BC_REFLECTION || 
		      BCtypeN[i] == BC_NO_SLIP ||
		      BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		      BCtypeN[i] == BC_MOVING_WALL ||
		      BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		      BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		      BCtypeN[i] == BC_BURNING_SURFACE ||
		      BCtypeN[i] == BC_MASS_INJECTION)) {
	    if (i > INl && i < INu) {
 	      norm_dir = HALF*(nfaceN(i, JCu) + 
                               nfaceN(i-1, JCu));
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		X_norm = ((Node[i][JNu].X - 
			   Node[i][JNu-GCell].X) * norm_dir) * norm_dir;
		X_tan = (Node[i][JNu].X - 
			 Node[i][JNu-GCell].X) - X_norm;
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 X_norm - X_tan );
	      }
	    } else if ( i == INl ) {
	      //norm_dir = nfaceN(i, JCu);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 (Node[i][JNu].X - 
					  Node[i][JNu-GCell].X) );
	      }
	    } else {
	      //norm_dir = nfaceN(i-1, JCu);
	      for(int GCell=1; GCell<=Nghost; ++GCell){
		Node[i][JNu+GCell].X = ( Node[i][JNu].X +
					 (Node[i][JNu].X - 
					  Node[i][JNu-GCell].X) );
	      }
	    } /* endif */
	  } else if (BCtypeN[i] == BCtypeN[i-1] &&
		     BCtypeN[i] == BC_PERIODIC) {
	    for(int GCell=1; GCell<=Nghost; ++GCell){
	      Node[i][JNu+GCell].X = ( Node[i][JNu].X +
				       (Node[i][JNl+GCell].X - 
					Node[i][JNl].X) );
	    }
	  } /* endif */
    } /* endfor */
    
    // Update SW and SE corner nodes.
    for (int j = JNl-Nghost ; j <= JNl ; ++j) {
      if (BCtypeW[j] != BC_REFLECTION &&
	  BCtypeW[j] != BC_PERIODIC &&  
	  BCtypeW[j] != BC_NO_SLIP &&
	  BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j] != BC_MOVING_WALL &&
	  BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j] != BC_BURNING_SURFACE &&
	  BCtypeW[j] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = Node[INl][j].X -
	    (Node[INl+GCell][j].X - 
	     Node[INl][j].X);
	}
      } else if (BCtypeW[j] == BC_REFLECTION || 
		 BCtypeW[j] == BC_NO_SLIP ||
		 BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeW[j] == BC_MOVING_WALL ||
		 BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeW[j] == BC_BURNING_SURFACE ||
		 BCtypeW[j] == BC_MASS_INJECTION) {
	if (j != JNl-Nghost) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	} else {
	  norm_dir = - nfaceW(ICl, j); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INl+GCell][j].X - 
		     Node[INl][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INl+GCell][j].X - 
		   Node[INl][j].X) - X_norm;
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   X_norm + X_tan );
	}
      } else if (BCtypeW[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } /* endif */
      
      if (BCtypeE[j] != BC_REFLECTION &&
	  BCtypeE[j] != BC_PERIODIC &&   
	  BCtypeE[j] != BC_NO_SLIP &&
	  BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j] != BC_MOVING_WALL &&
	  BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j] != BC_BURNING_SURFACE &&
	  BCtypeE[j] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j] == BC_REFLECTION ||
		 BCtypeE[j] == BC_NO_SLIP ||
		 BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeE[j] == BC_MOVING_WALL ||
		 BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeE[j] == BC_MOVING_WALL_HEATFLUX  ||
		 BCtypeE[j] == BC_BURNING_SURFACE ||
		 BCtypeE[j] == BC_MASS_INJECTION) {
	if (j != JNl-Nghost) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	} else {
	  norm_dir = nfaceE(ICu, j); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INu][j].X - 
		     Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INu][j].X - 
		   Node[INu-GCell][j].X) - X_norm;
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   X_norm - X_tan );
	}
      } else if (BCtypeE[j] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */
    
    // Update NW and NE corner nodes.
    for (int j = JNu ; j <= JNu+Nghost ; ++j) {
      if (BCtypeW[j-1] != BC_REFLECTION &&
	  BCtypeW[j-1] != BC_PERIODIC &&  
	  BCtypeW[j-1] != BC_NO_SLIP &&
	  BCtypeW[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeW[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeW[j-1] != BC_MOVING_WALL &&
	  BCtypeW[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeW[j-1] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeW[j-1] != BC_BURNING_SURFACE &&
	  BCtypeW[j-1] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } else if (BCtypeW[j-1] == BC_REFLECTION || 
		 BCtypeW[j-1] == BC_NO_SLIP ||
		 BCtypeW[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeW[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeW[j-1] == BC_MOVING_WALL ||
		 BCtypeW[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeW[j-1] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeW[j-1] == BC_BURNING_SURFACE ||
		 BCtypeW[j-1] == BC_MASS_INJECTION) {
	if (j != JNu+Nghost) {
	  norm_dir = - HALF*(nfaceW(ICl, j) + 
			     nfaceW(ICl, j-1));
	} else {
	  norm_dir = - nfaceW(ICl, j-1); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INl+GCell][j].X - 
		     Node[INl][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INl+GCell][j].X - 
		   Node[INl][j].X) - X_norm;
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   X_norm + X_tan );
	}
      } else if (BCtypeW[j-1] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INl-GCell][j].X = ( Node[INl][j].X -
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } /* endif */
      
      if (BCtypeE[j-1] != BC_REFLECTION &&
	  BCtypeE[j-1] != BC_PERIODIC && 
	  BCtypeE[j-1] != BC_NO_SLIP &&
	  BCtypeE[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	  BCtypeE[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
	  BCtypeE[j-1] != BC_MOVING_WALL &&
	  BCtypeE[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
	  BCtypeE[j-1] != BC_MOVING_WALL_HEATFLUX &&
	  BCtypeE[j-1] != BC_BURNING_SURFACE &&
	  BCtypeE[j-1] != BC_MASS_INJECTION) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INu][j].X - 
				    Node[INu-GCell][j].X) );
	}
      } else if (BCtypeE[j-1] == BC_REFLECTION || 
		 BCtypeE[j-1] == BC_NO_SLIP ||
		 BCtypeE[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		 BCtypeE[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		 BCtypeE[j-1] == BC_MOVING_WALL ||
		 BCtypeE[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		 BCtypeE[j-1] == BC_MOVING_WALL_HEATFLUX ||
		 BCtypeE[j-1] == BC_BURNING_SURFACE ||
		 BCtypeE[j-1] == BC_MASS_INJECTION) {
	if (j != JNu+Nghost) {
	  norm_dir = HALF*(nfaceE(ICu, j) + 
			   nfaceE(ICu, j-1));
	} else {
	  norm_dir = nfaceE(ICu, j-1); 
	} /* endif */
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  X_norm = ((Node[INu][j].X - 
		     Node[INu-GCell][j].X) * norm_dir) * norm_dir;
	  X_tan = (Node[INu][j].X - 
		   Node[INu-GCell][j].X) - X_norm;
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   X_norm - X_tan );
	}
      } else if (BCtypeE[j-1] == BC_PERIODIC) {
	for(int GCell=1; GCell<=Nghost; ++GCell){
	  Node[INu+GCell][j].X = ( Node[INu][j].X +
				   (Node[INl+GCell][j].X - 
				    Node[INl][j].X) );
	}
      } /* endif */
    } /* endfor */

    // Require update of the ghost cells geometric properties.
    Schedule_Ghost_Cells_Update();
    
    // Update corners
    Update_Corner_Ghost_Nodes();
}


/*!
 * Updates the corner ghost nodes for the quadrilateral 
 * mesh block.                                                                                                
 */
void Grid2D_Quad_Block_HO::Update_Corner_Ghost_Nodes(void) {
  
  Vector2D norm_dir, X_norm, X_tan;
  
  // SOUTH-WEST corner:
  if (BCtypeS[INl-1] == BC_NONE &&
      BCtypeW[JNl-1] == BC_NONE) {
    // Do nothing.
    
  } else if ((BCtypeS[INl-1] != BC_REFLECTION &&
	      BCtypeS[INl-1] != BC_PERIODIC &&
	      BCtypeS[INl-1] != BC_NO_SLIP &&
	      BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[INl-1] != BC_MOVING_WALL &&
	      BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	      BCtypeS[INl-1] != BC_MASS_INJECTION &&
	      BCtypeS[INl-1] != BC_RINGLEB_FLOW) &&
	     BCtypeW[JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	Node[i][JNl-ng].X = (Node[i][JNl].X +
			     (Node[i][JNl].X - 
			      Node[i][JNl+ng].X) );
      }
    }
    
  } else if ((BCtypeS[INl-1] == BC_NONE ||
	      (BCtypeS[INl-1] != BC_REFLECTION &&
	       BCtypeS[INl-1] != BC_PERIODIC &&
	       BCtypeS[INl-1] != BC_NO_SLIP &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeS[INl-1] != BC_MASS_INJECTION &&
	       BCtypeS[INl-1] != BC_RINGLEB_FLOW)) &&
	     (BCtypeW[JNl-1] != BC_REFLECTION &&
	      BCtypeW[JNl-1] != BC_PERIODIC &&  
	      BCtypeW[JNl-1] != BC_NO_SLIP &&
	      BCtypeW[JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeW[JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeW[JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeW[JNl-1] != BC_BURNING_SURFACE &&
	      BCtypeW[JNl-1] != BC_MASS_INJECTION)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
 	Node[INl-ng][j].X = ( Node[INl][j].X -
			      (Node[INl+ng][j].X - 
			       Node[INl][j].X) );
      }
    }
    
  } else if ((BCtypeS[INl-1] == BC_NONE ||
	      (BCtypeS[INl-1] != BC_REFLECTION &&
	       BCtypeS[INl-1] != BC_PERIODIC &&
	       BCtypeS[INl-1] != BC_NO_SLIP &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeS[INl-1] != BC_MASS_INJECTION &&
	       BCtypeS[INl-1] != BC_RINGLEB_FLOW)) &&
	     (BCtypeW[JNl-1] == BC_REFLECTION ||
	      BCtypeW[JNl-1] == BC_PERIODIC ||
	      BCtypeW[JNl-1] == BC_NO_SLIP ||
	      BCtypeW[JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeW[JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeW[JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeW[JNl-1] == BC_BURNING_SURFACE ||
	      BCtypeW[JNl-1] == BC_MASS_INJECTION)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
	if (j == JNl-Nghost) {
	  norm_dir = - nfaceW(ICl,JCl-Nghost);
	} else {
	  norm_dir = - HALF*(nfaceW(ICl,j) + 
			     nfaceW(ICl,j-1));
	}
	X_norm = ((Node[INl+ng][j].X - 
		   Node[INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INl+ng][j].X -
		 Node[INl][j].X) - X_norm;
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      X_norm + X_tan );
      }
    }

  } else if (BCtypeS[INl-1] == BC_REFLECTION ||
	     BCtypeS[INl-1] == BC_PERIODIC ||
	     BCtypeS[INl-1] == BC_NO_SLIP ||
	     BCtypeS[INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeS[INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeS[INl-1] == BC_MOVING_WALL ||
	     BCtypeS[INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeS[INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeS[INl-1] == BC_BURNING_SURFACE ||
	     BCtypeS[INl-1] == BC_MASS_INJECTION ||
	     BCtypeS[INl-1] == BC_RINGLEB_FLOW) {
    // Reflect cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      //for (int i = INl-Nghost; i <= INl; i++) {
      for (int i = INl; i >= INl-Nghost; i--) {
	if (i == INl-Nghost) {
	  norm_dir = - nfaceS(ICl-Nghost,JCl);
	} else {
	  norm_dir = - HALF*(nfaceS(i,JCl) +
			     nfaceS(i-1,JCl));
	}
	X_norm = ((Node[i][JNl+ng].X - 
		   Node[i][JNl].X)*norm_dir)*norm_dir;
	X_tan = (Node[i][JNl+ng].X - 
		 Node[i][JNl].X) - X_norm;
	Node[i][JNl-ng].X = ( Node[i][JNl].X -
			      X_norm + X_tan );
      }
    }
  }

  // SOUTH-EAST corner:
  if (BCtypeS[INu+1] == BC_NONE &&
      BCtypeE[JNl-1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeS[INu+1] != BC_REFLECTION &&
	      BCtypeS[INu+1] != BC_PERIODIC &&
	      BCtypeS[INu+1] != BC_NO_SLIP &&
	      BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeS[INu+1] != BC_MOVING_WALL &&
	      BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	      BCtypeS[INu+1] != BC_MASS_INJECTION &&
	      BCtypeS[INu+1] != BC_RINGLEB_FLOW) &&
	     BCtypeE[JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	Node[i][JNl-ng].X = ( Node[i][JNl].X +
			      (Node[i][JNl].X - 
			       Node[i][JNl+ng].X) );
      }
    }

  } else if ((BCtypeS[INu-1] == BC_NONE ||
	      (BCtypeS[INu+1] != BC_REFLECTION &&
	       BCtypeS[INu+1] != BC_PERIODIC &&
	       BCtypeS[INu+1] != BC_NO_SLIP &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INu+1] != BC_MOVING_WALL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeS[INu+1] != BC_MASS_INJECTION &&
	       BCtypeS[INu+1] != BC_RINGLEB_FLOW)) &&
	     (BCtypeE[JNl-1] != BC_REFLECTION &&
	      BCtypeE[JNl-1] != BC_PERIODIC &&  
	      BCtypeE[JNl-1] != BC_NO_SLIP &&
	      BCtypeE[JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeE[JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeE[JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeE[JNl-1] != BC_BURNING_SURFACE &&
	      BCtypeE[JNl-1] != BC_MASS_INJECTION)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl-Nghost; j <= JNl; j++) {
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      (Node[INu][j].X - 
			       Node[INu-ng][j].X) );
      }
    }

  } else if ((BCtypeS[INu+1] == BC_NONE ||
	      (BCtypeS[INu+1] != BC_REFLECTION &&
	       BCtypeS[INu+1] != BC_PERIODIC &&
	       BCtypeS[INu+1] != BC_NO_SLIP &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeS[INu+1] != BC_MOVING_WALL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeS[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeS[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeS[INu+1] != BC_MASS_INJECTION &&
	       BCtypeS[INu+1] != BC_RINGLEB_FLOW)) &&
	     (BCtypeE[JNl-1] == BC_REFLECTION ||
	      BCtypeE[JNl-1] == BC_PERIODIC ||
	      BCtypeE[JNl-1] == BC_NO_SLIP ||
	      BCtypeE[JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeE[JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeE[JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeE[JNl-1] == BC_BURNING_SURFACE ||
	      BCtypeE[JNl-1] == BC_MASS_INJECTION)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNl; j >= JNl-Nghost; j--) {
	if (j == JNl-Nghost) {
	  norm_dir = nfaceE(ICu,JCl-Nghost);
	} else {
	  norm_dir = HALF*(nfaceE(ICu,j) + 
			   nfaceE(ICu,j+1));
	}
	norm_dir = nfaceE(ICu,j);
	X_norm = ((Node[INu][j].X - 
		   Node[INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INu][j].X - 
		 Node[INu-ng][j].X) - X_norm;
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      X_norm - X_tan );
      }
    }

  } else if (BCtypeS[INu+1] == BC_REFLECTION ||
	     BCtypeS[INu+1] == BC_PERIODIC ||
	     BCtypeS[INu+1] == BC_NO_SLIP ||
	     BCtypeS[INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeS[INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeS[INu+1] == BC_MOVING_WALL ||
	     BCtypeS[INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeS[INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeS[INu+1] == BC_BURNING_SURFACE ||
	     BCtypeS[INu+1] == BC_MASS_INJECTION ||
	     BCtypeS[INu+1] == BC_RINGLEB_FLOW) {
    // Reflect cells south.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	if (i == INu+Nghost) {
	  norm_dir = -nfaceS(ICu+Nghost,JCl);
	} else {
	  norm_dir = -HALF*(nfaceS(i,JCl) +
			    nfaceS(i-1,JCl));
	}
 	X_norm = ((Node[i][JNl+ng].X - 
 		   Node[i][JNl].X)*norm_dir)*norm_dir;
 	X_tan = (Node[i][JNl+ng].X - 
 		 Node[i][JNl].X) - X_norm;
  	Node[i][JNl-ng].X = ( Node[i][JNl].X -
			      X_norm + X_tan );
      }
    }
  }

  // NORTH-WEST corner:
  if (BCtypeN[INl-1] == BC_NONE &&
      BCtypeW[JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeN[INl-1] != BC_REFLECTION &&
	      BCtypeN[INl-1] != BC_PERIODIC &&
	      BCtypeN[INl-1] != BC_NO_SLIP &&
	      BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[INl-1] != BC_MOVING_WALL &&
	      BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	      BCtypeN[INl-1] != BC_MASS_INJECTION) &&
	     BCtypeW[JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      (Node[i][JNu].X - 
			       Node[i][JNu-ng].X) );
      }
    }

  } else if ((BCtypeN[INl-1] == BC_NONE ||
	      (BCtypeN[INl-1] != BC_REFLECTION &&
	       BCtypeN[INl-1] != BC_PERIODIC &&
	       BCtypeN[INl-1] != BC_NO_SLIP &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INl-1] != BC_MOVING_WALL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeN[INl-1] != BC_MASS_INJECTION)) &&
	     (BCtypeW[JNu+1] != BC_REFLECTION &&
	      BCtypeW[JNu+1] != BC_PERIODIC &&  
	      BCtypeW[JNu+1] != BC_NO_SLIP &&
	      BCtypeW[JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeW[JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeW[JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeW[JNu+1] != BC_BURNING_SURFACE &&
	      BCtypeW[JNu+1] != BC_MASS_INJECTION)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      (Node[INl+ng][j].X - 
			       Node[INl][j].X) );
      }
    }

  } else if ((BCtypeN[INl-1] == BC_NONE ||
	      (BCtypeN[INl-1] != BC_REFLECTION &&
	       BCtypeN[INl-1] != BC_PERIODIC &&
	       BCtypeN[INl-1] != BC_NO_SLIP &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INl-1] != BC_MOVING_WALL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INl-1] != BC_BURNING_SURFACE &&
	       BCtypeN[INl-1] != BC_MASS_INJECTION)) &&
	     (BCtypeW[JNu+1] == BC_REFLECTION ||
	      BCtypeW[JNu+1] == BC_PERIODIC ||
	      BCtypeW[JNu+1] == BC_NO_SLIP ||
	      BCtypeW[JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeW[JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeW[JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeW[JNu+1] == BC_BURNING_SURFACE ||
	      BCtypeW[JNu+1] == BC_MASS_INJECTION)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	if (j == JNu+Nghost) {
	  norm_dir = - nfaceW(ICl,JCu+Nghost);
	} else {
	  norm_dir = - HALF*(nfaceW(ICl,j) + 
			     nfaceW(ICl,j-1));
	}
	X_norm = ((Node[INl+ng][j].X -
		   Node[INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INl+ng][j].X - 
		 Node[INl][j].X) - X_norm;
	Node[INl-ng][j].X = ( Node[INl][j].X -
			      X_norm + X_tan );
      }
    }

  } else if (BCtypeN[INl-1] == BC_REFLECTION ||
	     BCtypeN[INl-1] == BC_PERIODIC ||
	     BCtypeN[INl-1] == BC_NO_SLIP ||
	     BCtypeN[INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeN[INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeN[INl-1] == BC_MOVING_WALL ||
	     BCtypeN[INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeN[INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeN[INl-1] == BC_BURNING_SURFACE ||
	     BCtypeN[INl-1] == BC_MASS_INJECTION) {
    // Reflect cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INl-Nghost; i <= INl; i++) {
	if (i == INl-Nghost) {
	  norm_dir = -nfaceN(ICl-Nghost,JCu);
	} else {
	  norm_dir = -HALF*(nfaceN(i,JCu) + 
			    nfaceN(i-1,JCu));
	}
	X_norm = ((Node[i][JNu].X - 
		   Node[i][JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Node[i][JNu].X - 
		 Node[i][JNu-ng].X) - X_norm;
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      X_norm - X_tan );
      }
    }
  }

  // NORTH-EAST corner:
  if (BCtypeN[INu+1] == BC_NONE &&
      BCtypeE[JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((BCtypeN[INu+1] != BC_REFLECTION &&
	      BCtypeN[INu+1] != BC_PERIODIC &&
	      BCtypeN[INu+1] != BC_NO_SLIP &&
	      BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeN[INu+1] != BC_MOVING_WALL &&
	      BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	      BCtypeN[INu+1] != BC_MASS_INJECTION) &&
	     BCtypeE[JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      (Node[i][JNu].X - 
			       Node[i][JNu-ng].X) );
      }
    }

  } else if ((BCtypeN[INu+1] == BC_NONE ||
	      (BCtypeN[INu+1] != BC_REFLECTION &&
	       BCtypeN[INu+1] != BC_PERIODIC &&
	       BCtypeN[INu+1] != BC_NO_SLIP &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INu+1] != BC_MOVING_WALL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeN[INu+1] != BC_MASS_INJECTION)) &&
	     (BCtypeE[JNu+1] != BC_REFLECTION &&
	      BCtypeE[JNu+1] != BC_PERIODIC &&  
	      BCtypeE[JNu+1] != BC_NO_SLIP &&
	      BCtypeE[JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      BCtypeE[JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      BCtypeE[JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      BCtypeE[JNu+1] != BC_BURNING_SURFACE &&
	      BCtypeE[JNu+1] != BC_MASS_INJECTION)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      (Node[INu][j].X - 
			       Node[INu-ng][j].X) );
      }
    }

  } else if ((BCtypeN[INu+1] == BC_NONE ||
	      (BCtypeN[INu+1] != BC_REFLECTION &&
	       BCtypeN[INu+1] != BC_PERIODIC &&
	       BCtypeN[INu+1] != BC_NO_SLIP &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       BCtypeN[INu+1] != BC_MOVING_WALL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       BCtypeN[INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       BCtypeN[INu+1] != BC_BURNING_SURFACE &&
	       BCtypeN[INu+1] != BC_MASS_INJECTION)) &&
	     (BCtypeE[JNu+1] == BC_REFLECTION ||
	      BCtypeE[JNu+1] == BC_PERIODIC ||
	      BCtypeE[JNu+1] == BC_NO_SLIP ||
	      BCtypeE[JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      BCtypeE[JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      BCtypeE[JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      BCtypeE[JNu+1] == BC_BURNING_SURFACE ||
	      BCtypeE[JNu+1] == BC_MASS_INJECTION)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int j = JNu; j <= JNu+Nghost; j++) {
	if (j == JNu+Nghost) {
	  norm_dir = nfaceE(ICu,JCu+Nghost);
	} else {
	  norm_dir = HALF*(nfaceE(ICu,j) + 
			   nfaceE(ICu,j-1));
	}
	X_norm = ((Node[INu][j].X - 
		   Node[INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Node[INu][j].X - 
		 Node[INu-ng][j].X) - X_norm;
	Node[INu+ng][j].X = ( Node[INu][j].X +
			      X_norm - X_tan );
      }
    }

  } else if (BCtypeN[INu+1] == BC_REFLECTION ||
	     BCtypeN[INu+1] == BC_PERIODIC ||
	     BCtypeN[INu+1] == BC_NO_SLIP ||
	     BCtypeN[INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     BCtypeN[INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     BCtypeN[INu+1] == BC_MOVING_WALL ||
	     BCtypeN[INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     BCtypeN[INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     BCtypeN[INu+1] == BC_BURNING_SURFACE ||
	     BCtypeN[INu+1] == BC_MASS_INJECTION) {
    // Reflect cells north.
    for (int ng = 1; ng <= Nghost; ng++) {
      for (int i = INu; i <= INu+Nghost; i++) {
	if (i == INu+Nghost) {
	  norm_dir = -nfaceN(ICu+Nghost,JCu);
	} else {
	  norm_dir = -HALF*(nfaceN(i,JCu) + 
			    nfaceN(i-1,JCu));
	}
	X_norm = ((Node[i][JNu].X - 
		   Node[i][JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Node[i][JNu].X - 
		 Node[i][JNu-ng].X) - X_norm;
	Node[i][JNu+ng].X = ( Node[i][JNu].X +
			      X_norm - X_tan );
      }
    }
  }

  // Require the update of the corner ghost cells' geometric properties
  Schedule_Corner_Ghost_Cells_Update();

  // Update cell info.
  Update_Cells();
}


/*!
 * Updates the cell information for the quadrilateral
 * mesh block. Different levels of update are possible
 * based on the corresponding flag values.
 * If all the flags are OFF this subroutine does nothing.
 */
void Grid2D_Quad_Block_HO::Update_Cells(void) {

  // Decide which level of update is required
  if (InteriorMeshUpdate == ON && GhostCellsUpdate == ON ){
    // Update all the cells and the boundary spline information
    Update_All_Cells();
  } else if (InteriorMeshUpdate == ON) {
    // Update only the information of the interior cells and the boundary spline information
    Update_Interior_Cells();
  } else if (GhostCellsUpdate == ON) {
    // Update only the information of the ghost cells
    Update_Ghost_Cells();
  } else if (CornerGhostCellsUpdate == ON){
    // Update only the information of the corner ghost cells
    Update_Corner_Ghost_Cells();
  }
}

/*!
 * Updates the cell information for the quadrilateral mesh block.
 * (i.e. all cells).
 */
void Grid2D_Quad_Block_HO::Update_All_Cells(void) {

  int i,j;

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
      Update_Cell(i,j);
    } /* endfor */
  } /* endfor */

  if ( HighOrderBoundaryRepresentation == OFF ){
    // Confirm the update
    Confirm_Mesh_Update_Everywhere();
    return;
  }

}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * interior cells (i.e. no ghost cells).
 */
void Grid2D_Quad_Block_HO::Update_Interior_Cells(void) {

  int i,j;

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl; j <= JCu ; ++j) {
    for ( i = ICl; i <= ICu; ++i) {
      Update_Cell(i,j);
    } /* endfor */
  } /* endfor */
  
  if ( HighOrderBoundaryRepresentation == OFF ){
    // Confirm the update
    Confirm_Mesh_Update_Everywhere();
    return;
  }

}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * ghost cells (i.e. no interior cells and no boundary splines).
 */
void Grid2D_Quad_Block_HO::Update_Ghost_Cells(void) {

  int i,j;

  // Update cell information assuming straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
    // Update the West face ghost cells
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the East face ghost cells
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  for (i = ICl; i <= ICu; ++i){
    // Update the South face ghost cells
    for (j = JCl-Nghost; j<=JCl-1; ++j){
      Update_Cell(i,j);
    } // endif(j)
    
    // Update the North face ghost cells
    for (j = JCu+1; j<=JCu+Nghost; ++j){
      Update_Cell(i,j);
    } // endif(j)
  } // endif(i)
  
  if ( HighOrderBoundaryRepresentation == OFF ){
    // Confirm the update
    Confirm_Mesh_Update_Everywhere();
    return;
  }

  if ( (BndNorthSpline.Xp == NULL || BndNorthSpline.bc[0] == BC_NONE || BndNorthSpline.bc[0] == BC_PERIODIC) &&
       (BndSouthSpline.Xp == NULL || BndSouthSpline.bc[0] == BC_NONE || BndSouthSpline.bc[0] == BC_PERIODIC) &&
       (BndEastSpline.Xp == NULL  || BndEastSpline.bc[0] == BC_NONE || BndEastSpline.bc[0] == BC_PERIODIC) && 
       (BndWestSpline.Xp == NULL  || BndWestSpline.bc[0] == BC_NONE || BndWestSpline.bc[0] == BC_PERIODIC) ){   
    // No curved boundaries

    // Confirm the update
    Confirm_Mesh_Update_Everywhere();
    return;
  }

  /* Recompute the geometric properties of those interior cells that are near a curved boundary.
     Obs1. The "SPLINE2D_LINEAR" is also considered a curved boundary because it might have a 
     sharp point between two nodes.
  */

}

/*!
 * Updates the cell information for the quadrilateral mesh block 
 * corner ghost cells (i.e. no interior cells, no ghost cells other
 * than the corner ones, no boundary splines).
 */
void Grid2D_Quad_Block_HO::Update_Corner_Ghost_Cells(void) {

  int i,j;

  // Update cell information with straight boundaries (i.e. every edge of the cell is a line segment)
  for ( j = JCl-Nghost ; j <= JCl-1 ; ++j) {
    // Update the South-West corner
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the South-East corner
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  for ( j = JCu+1 ; j <= JCu+Nghost ; ++j) {
    // Update the North-West corner
    for (i = ICl-Nghost; i <= ICl-1; ++i){
      Update_Cell(i,j);
    } // endif (i)

    // Update the North-East corner
    for (i = ICu+1; i <= ICu+Nghost; ++i){
      Update_Cell(i,j);
    } // endif (i)
  } // endif(j)

  // Confirm the update
  Confirm_Ghost_Cells_Update();
}

/*!
 * Check the validity of the quadrilateral mesh block.  
 * Returns a non-zero result if mesh is not valid.      
 *                                                      
 */
int Grid2D_Quad_Block_HO::Check_Quad_Block(void) {

    int i, j;

    for ( j = JCl ; j <= JCu ; ++j) {
        for ( i = ICl ; i <= ICu ; ++i) {
	  if (Cell[i][j].A <= ZERO) {
	    cout << endl << i << " " << j;
	    cout << endl << ICl << " " << ICu;
	    cout << endl << JCl << " " << JCu;
	    cout << endl << Cell[i][j].A;
	    cout << endl << Cell[i][j].Xc;
	    cout << endl << Node[i][j].X;
	    cout << endl << Node[i+1][j].X;
	    cout << endl << Node[i][j+1].X;
	    cout << endl << Node[i+1][j+1].X;
	    cout.flush();
   	    return(1);
          } /* endif */
        } /* endfor */
    } /* endfor */

    return(0);

}

/*!
 * Writes definition file information for 2D            
 * quadrilateral grid block to the specified output     
 * stream for retrieval and re-use purposes.            
 *                                                      
 */
void Grid2D_Quad_Block_HO::Write_Quad_Block_Definition(ostream &Out_File) {

    Out_File << setprecision(14);
    if (NNi == 0 || NNj == 0) {
       Out_File << NCi << " " 
	        << NCj << " "
		<< Nghost << "\n";
    } else {
       Out_File << NCi << " " 
	        << NCj << " "
		<< Nghost << "\n" 
                << BndNorthSpline.np << "\n"
                << BndNorthSpline.type << "\n"
                << BndNorthSpline 
                << BndSouthSpline.np << "\n"
                << BndSouthSpline.type << "\n"
                << BndSouthSpline
                << BndEastSpline.np << "\n"
                << BndEastSpline.type << "\n"
                << BndEastSpline 
                << BndWestSpline.np << "\n"
                << BndWestSpline.type << "\n"
                << BndWestSpline 
                << GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH << "\n";
       Out_File.setf(ios::scientific);
       Out_File << StretchI << " " 
                << BetaI << " " 
                << TauI << "\n"
                << StretchJ << " " 
                << BetaJ << " " 
                << TauJ << "\n";
       Out_File.unsetf(ios::scientific);
       Out_File << OrthogonalN << " "
                << OrthogonalS << " "
                << OrthogonalE << " "
                << OrthogonalW << "\n";
    } /* endif */
    Out_File << setprecision(6);

}


/*!
 * Reads definition file for 2D quadrilateral grid      
 * block from the specified output stream.              
 *                                                      
 */
void Grid2D_Quad_Block_HO::Read_Quad_Block_Definition(istream &In_File) {
  
  int i, j, ng, k, kx, ky, kx_max, ky_max, npts, spline_type;
  int node_init_procedure;
  double S_i, S_j, 
    s_north, s_south, s_east, s_west,
    smax_north, smax_south, smax_east, smax_west, 
    s_i, s_j, smax_i, smax_j,
    w_north, w_south, w_east, w_west, w_total;
  double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
    S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
  Vector2D x_north, x_south, x_east, x_west;
  Spline2D_HO S1_NS, S2_NS, S1_EW, S2_EW;

  In_File.setf(ios::skipws);
  In_File >> i >> j >> ng;
  In_File.unsetf(ios::skipws);

  if (i != 0 && j != 0 && ng != 0) {
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block as required. */
    allocate(i, j, ng);
    
    /* For each of the north, south, east, and west boundaries
       of this mesh block, read in the number of spline points, 
       allocate memory for the boundary splines, read in the 
       spline points, and finally calculate the spline 
       pathlengths. */

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    if (BndNorthSpline.np != 0) BndNorthSpline.deallocate();
    BndNorthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndNorthSpline.settype(spline_type);
    In_File >> BndNorthSpline;
    BndNorthSpline.pathlength();
    SminN = BndNorthSpline.sp[0];
    SmaxN = BndNorthSpline.sp[BndNorthSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    if (BndSouthSpline.np != 0) BndSouthSpline.deallocate();
    BndSouthSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndSouthSpline.settype(spline_type);
    In_File >> BndSouthSpline;
    BndSouthSpline.pathlength();
    SminS = BndSouthSpline.sp[0];
    SmaxS = BndSouthSpline.sp[BndSouthSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    if (BndEastSpline.np != 0) BndEastSpline.deallocate();
    BndEastSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndEastSpline.settype(spline_type);
    In_File >> BndEastSpline;
    BndEastSpline.pathlength();
    SminE = BndEastSpline.sp[0];
    SmaxE = BndEastSpline.sp[BndEastSpline.np-1];

    In_File.setf(ios::skipws);
    In_File >> npts;
    In_File.unsetf(ios::skipws);
    if (BndWestSpline.np != 0) BndWestSpline.deallocate();
    BndWestSpline.allocate(npts);
    In_File.setf(ios::skipws);
    In_File >> spline_type;
    In_File.unsetf(ios::skipws);
    BndWestSpline.settype(spline_type);
    In_File >> BndWestSpline;
    BndWestSpline.pathlength();
    SminW = BndWestSpline.sp[0];
    SmaxW = BndWestSpline.sp[BndWestSpline.np-1];

    /* Read the node initialization procedure for this 
       quadrilateral grid block. */

    In_File.setf(ios::skipws);
    In_File >> node_init_procedure;
    In_File.unsetf(ios::skipws);
      
    /* Input the type of node distribution stretching 
       functions to be used for each of the coordinate
       directions. Also read in the related stretching
       parameters. */

    In_File.setf(ios::skipws);
    In_File >> StretchI >> BetaI >> TauI;
    In_File >> StretchJ >> BetaJ >> TauJ;
    In_File.unsetf(ios::skipws);

    /* Input the grid orthogonality specifiers for
       the north, south, east, and west boundaries. */

    In_File.setf(ios::skipws);
    In_File >> OrthogonalN >> OrthogonalS
	    >> OrthogonalE >> OrthogonalW;
    In_File.unsetf(ios::skipws);

    /* Compute the interior nodes for the quadrilateral mesh block. */

    smax_north = BndNorthSpline.sp[BndNorthSpline.np-1]-
      BndNorthSpline.sp[0];
    smax_south = BndSouthSpline.sp[BndSouthSpline.np-1]-
      BndSouthSpline.sp[0];

    smax_east = BndEastSpline.sp[BndEastSpline.np-1]-
      BndEastSpline.sp[0];
    smax_west = BndWestSpline.sp[BndWestSpline.np-1]-
      BndWestSpline.sp[0];

    dS_i = ONE/(TEN*double(INu-INl));
    dS_j = ONE/(TEN*double(JNu-JNl));

    for ( j = JNl ; j <= JNu ; ++j) {
      for ( i = INl ; i <= INu ; ++i) {
	S_j = double(j-JNl)/double(JNu-JNl);
	S_j = StretchingFcn(S_j, BetaJ, TauJ, StretchJ);
	s_east = S_j * smax_east + BndEastSpline.sp[0];
	x_east = Spline(s_east, BndEastSpline);
	s_west = S_j * smax_west + BndWestSpline.sp[0];
	x_west = Spline(s_west, BndWestSpline);

	S_i = double(i-INl)/double(INu-INl);
	S_i = StretchingFcn(S_i, BetaI, TauI, StretchI);
	s_north = S_i * smax_north + BndNorthSpline.sp[0];
	x_north = Spline(s_north, BndNorthSpline);
	s_south = S_i * smax_south + BndSouthSpline.sp[0];
	x_south = Spline(s_south, BndSouthSpline);

	if (j == JNl) {
	  Node[i][j].X = x_south;
	} else if (j == JNu) {
	  Node[i][j].X = x_north;
	} else if (i == INl) {
	  Node[i][j].X = x_west;
	} else if (i == INu) {
	  Node[i][j].X = x_east;
	} else {
	  s_i = (ONE-S_j)*s_south + S_j*s_north;
	  s_j = (ONE-S_i)*s_west + S_i*s_east;

	  smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
	  smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

	  switch(node_init_procedure) {
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
	    w_south = (ONE-s_j/smax_j);
	    w_north = (s_j/smax_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
	    w_west =  (ONE-s_i/smax_i);
	    w_east =  (s_i/smax_i);
	    w_total = w_east + w_west;
	    Node[i][j].X = (w_west*x_west + 
				 w_east*x_east)/w_total;
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
	    kx_max = 5;
	    dS_i = ONE/double(kx_max-1);
	    kx = int(floor(S_i/dS_i));
	    S1_i = max(ZERO, double(kx)*dS_i);
	    S2_i = min(ONE, S1_i + dS_i);

	    ky_max = 7;
	    dS_j = ONE/double(ky_max-1);
	    ky = int(floor(S_j/dS_j));
	    S1_j = max(ZERO, double(ky)*dS_j);
	    S2_j = min(ONE, S1_j + dS_j);

	    if (S1_i - TOLER > ZERO) {
	      S1_NS.allocate(2);
	      S1_NS.settype(SPLINE2D_LINEAR);

	      s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	      S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S1_i * smax_north + BndNorthSpline.sp[0],
		S1_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	      S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S1_NS.pathlength();
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } else {
	      S1_NS = BndWestSpline;
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } /* endif */

	    if (S2_i + TOLER < ONE) {
	      S2_NS.allocate(2);
	      S2_NS.settype(SPLINE2D_LINEAR);

	      s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	      S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S2_i * smax_north + BndNorthSpline.sp[0],
		S2_NS.Xp[1] = Spline(s_north, BndNorthSpline);
	      S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S2_NS.pathlength();
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } else {
	      S2_NS = BndEastSpline;
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } /* endif */

	    if (S1_j - TOLER > ZERO) {
	      S1_EW.allocate(kx_max);
	      S1_EW.settype(SPLINE2D_LINEAR);

	      s_west = S1_j * smax_west + BndWestSpline.sp[0];
	      S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S1_j * smax_east + BndEastSpline.sp[0];
	      S1_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	      S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
		s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
		x_north = Spline(s_north, BndNorthSpline);

		s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
		x_south = Spline(s_south, BndSouthSpline);

		w_south =  ONE-S1_j;
		w_north =  S1_j;
		w_total = w_north + w_south;
		S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
		S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S1_EW.pathlength();
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } else {
	      S1_EW = BndSouthSpline;
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } /* endif */

	    if (S2_j + TOLER < ONE) {
	      S2_EW.allocate(kx_max);
	      S2_EW.settype(SPLINE2D_LINEAR);

	      s_west = S2_j * smax_west + BndWestSpline.sp[0];
	      S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S2_j * smax_east + BndEastSpline.sp[0];
	      S2_EW.Xp[kx_max-1] = Spline(s_east, BndEastSpline);
	      S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
		s_north = double(k)*smax_north/double(kx_max-1) + BndNorthSpline.sp[0];
		x_north = Spline(s_north, BndNorthSpline);

		s_south = double(k)*smax_south/double(kx_max-1) + BndSouthSpline.sp[0];
		x_south = Spline(s_south, BndSouthSpline);

		w_south =  ONE-S2_j;
		w_north =  S2_j;
		w_total = w_north + w_south;
		S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
		S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S2_EW.pathlength();
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } else {
	      S2_EW = BndNorthSpline;
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } /* endif */

	    s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	    x_north = Spline(s_north, S2_EW);
	    s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	    x_south = Spline(s_south, S1_EW);

	    s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	    x_east = Spline(s_east, S2_NS);
	    s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	    x_west = Spline(s_west, S1_NS);
      
	    if (ky == 0 || ky == ky_max - 2) {
	      w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	      w_north =  (S_j-S1_j)/(S2_j-S1_j);
	      w_total = w_south + w_north;
	      Node[i][j].X = (w_south*x_south + 
				   w_north*x_north)/w_total;
	    } else {
	      w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	      w_east =  (S_i-S1_i)/(S2_i-S1_i);
	      w_total = w_east + w_west;
	      Node[i][j].X = (w_west*x_west + 
				   w_east*x_east)/w_total;
	    } /* endif */

	    S1_NS.deallocate(); S2_NS.deallocate();
	    S2_EW.deallocate(); S2_EW.deallocate();
	    break;
	    //============================================================
	  case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
	    kx_max = 7;
	    dS_i = ONE/double(kx_max-1);
	    kx = int(floor(S_i/dS_i));
	    S1_i = max(ZERO, double(kx)*dS_i);
	    S2_i = min(ONE, S1_i + dS_i);

	    ky_max = 5;
	    dS_j = ONE/double(ky_max-1);
	    ky = int(floor(S_j/dS_j));
	    S1_j = max(ZERO, double(ky)*dS_j);
	    S2_j = min(ONE, S1_j + dS_j);

	    if (S1_i - TOLER > ZERO) {
	      S1_NS.allocate(ky_max);
	      S1_NS.settype(SPLINE2D_LINEAR);

	      s_south = S1_i * smax_south + BndSouthSpline.sp[0];
	      S1_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S1_i * smax_north + BndNorthSpline.sp[0],
		S1_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	      S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
		s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
		x_east = Spline(s_east, BndEastSpline);

		s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
		x_west = Spline(s_west, BndWestSpline);

		w_west =  ONE-S1_i;
		w_east =  S1_i;
		w_total = w_east + w_west;
		S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
		S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S1_NS.pathlength();
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } else {
	      S1_NS = BndWestSpline;
	      smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
		S1_NS.sp[0];
	    } /* endif */

	    if (S2_i + TOLER < ONE) {
	      S2_NS.allocate(ky_max);
	      S2_NS.settype(SPLINE2D_LINEAR);

	      s_south = S2_i * smax_south + BndSouthSpline.sp[0];
	      S2_NS.Xp[0] = Spline(s_south, BndSouthSpline);
	      S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_north = S2_i * smax_north + BndNorthSpline.sp[0],
		S2_NS.Xp[ky_max-1] = Spline(s_north, BndNorthSpline);
	      S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

	      for ( k = 1 ; k <= ky_max-2 ; ++k) {
		s_east = double(k)*smax_east/double(ky_max-1) + BndEastSpline.sp[0];
		x_east = Spline(s_east, BndEastSpline);

		s_west = double(k)*smax_west/double(ky_max-1) + BndWestSpline.sp[0];
		x_west = Spline(s_west, BndWestSpline);

		w_west =  ONE-S2_i;
		w_east =  S2_i;
		w_total = w_east + w_west;
		S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
		S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
	      } /* endfor */

	      S2_NS.pathlength();
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } else {
	      S2_NS = BndEastSpline;
	      smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
		S2_NS.sp[0];
	    } /* endif */

	    if (S1_j - TOLER > ZERO) {
	      S1_EW.allocate(2);
	      S1_EW.settype(SPLINE2D_LINEAR);

	      s_west = S1_j * smax_west + BndWestSpline.sp[0];
	      S1_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S1_j * smax_east + BndEastSpline.sp[0];
	      S1_EW.Xp[1] = Spline(s_east, BndEastSpline);
	      S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S1_EW.pathlength();
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } else {
	      S1_EW = BndSouthSpline;
	      smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
		S1_EW.sp[0];
	    } /* endif */

	    if (S2_j + TOLER < ONE) {
	      S2_EW.allocate(2);
	      S2_EW.settype(SPLINE2D_LINEAR);
   
	      s_west = S2_j * smax_west + BndWestSpline.sp[0];
	      S2_EW.Xp[0] = Spline(s_west, BndWestSpline);
	      S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

	      s_east = S2_j * smax_east + BndEastSpline.sp[0];
	      S2_EW.Xp[1] = Spline(s_east, BndEastSpline);
	      S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

	      S2_EW.pathlength();
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } else {
	      S2_EW = BndNorthSpline;
	      smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
		S2_EW.sp[0];
	    } /* endif */

	    s_north = S_i * smax_S2_EW + S2_EW.sp[0];
	    x_north = Spline(s_north, S2_EW);
	    s_south = S_i * smax_S1_EW + S1_EW.sp[0];
	    x_south = Spline(s_south, S1_EW);

	    s_east = S_j * smax_S2_NS + S2_NS.sp[0];
	    x_east = Spline(s_east, S2_NS);
	    s_west = S_j * smax_S1_NS + S1_NS.sp[0];
	    x_west = Spline(s_west, S1_NS);
      
	    if (kx == 0 || kx == kx_max - 2) {
	      w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
	      w_east =  (S_i-S1_i)/(S2_i-S1_i);
	      w_total = w_east + w_west;
	      Node[i][j].X = (w_west*x_west + 
				   w_east*x_east)/w_total;
	    } else {
	      w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
	      w_north =  (S_j-S1_j)/(S2_j-S1_j);
	      w_total = w_south + w_north;
	      Node[i][j].X = (w_south*x_south + 
				   w_north*x_north)/w_total;
	    } /* endif */

	    S1_NS.deallocate(); S2_NS.deallocate();
	    S2_EW.deallocate(); S2_EW.deallocate();
	    break;
	    //============================================================
	  default:
	    w_south = (ONE-s_j/smax_j);
	    w_north = (s_j/smax_j);
	    w_total = w_south + w_north;
	    Node[i][j].X = (w_south*x_south + 
				 w_north*x_north)/w_total;
	    break;
	  } /* endswitch */

	} /* endif */

      } /* endfor */
    }/* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types at the quadrilateral 
       grid block boundaries. */

    Set_BCs();

    /* Compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for the quadrilateral mesh block. */

    Update_Cells();
  } /* endif */

}


/*!
 * Writes 2D quadrilateral grid block to the            
 * specified output stream for retrieval and re-use     
 * purposes.                                            
 */
void Grid2D_Quad_Block_HO::Write_Quad_Block(ostream &Out_File) {

  Out_File << setprecision(14) << *this << setprecision(6);

}

/*!
 * Reads 2D quadrilateral grid block from the           
 * specified input stream.                              
 *                                                      
 */
void Grid2D_Quad_Block_HO::Read_Quad_Block(istream &In_File) {

    In_File >> *this;

}

/*!
 * Translates or shifts the positions of the nodes of a 
 * quadrilateral grid block.                            
 *                                                      
 */
void Grid2D_Quad_Block_HO::Translate_Quad_Block(const Vector2D &V) {

    int i, j;

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           Node[i][j].X += V;
       } /* endfor */
    } /* endfor */

    for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
        for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
  	    Cell[i][j].Xc = centroid(i, j);
            Cell[i][j].A = area(i, j);
        } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Translate_Spline(V);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Translate_Spline(V);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Translate_Spline(V);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Translate_Spline(V);
 
}


/*!
 * Scales the quadrilateral grid block.                 
 *                                                      
 */
void Grid2D_Quad_Block_HO::Scale_Quad_Block(const double &Scaling_Factor) {

    int i, j;

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           Node[i][j].X = Node[i][j].X*Scaling_Factor;
       } /* endfor */
    } /* endfor */

    for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
        for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
  	    Cell[i][j].Xc = centroid(i, j);
            Cell[i][j].A = area(i, j);
        } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Scale_Spline(Scaling_Factor);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Scale_Spline(Scaling_Factor);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Scale_Spline(Scaling_Factor);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Scale_Spline(Scaling_Factor);

    SminN = SminN*Scaling_Factor;
    SmaxN = SmaxN*Scaling_Factor;
    SminS = SminS*Scaling_Factor;
    SmaxS = SmaxS*Scaling_Factor;
    SminE = SminE*Scaling_Factor;
    SmaxE = SmaxE*Scaling_Factor;
    SminW = SminW*Scaling_Factor;
    SmaxW = SmaxW*Scaling_Factor;
 
}

 
/*!
 * Rotates the quadrilateral grid block.                
 *                                                      
 */
void Grid2D_Quad_Block_HO::Rotate_Quad_Block(const double &Angle) {

    int i, j;
    double cos_angle, sin_angle;
    Vector2D X;

    cos_angle = cos(-Angle); 
    sin_angle = sin(-Angle);

    for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
       for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
           X.x = Node[i][j].X.x*cos_angle +
                 Node[i][j].X.y*sin_angle;
           X.y = - Node[i][j].X.x*sin_angle +
                   Node[i][j].X.y*cos_angle;
           Node[i][j].X = X;
       } /* endfor */
    } /* endfor */

    for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
        for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
  	    Cell[i][j].Xc = centroid(i, j);
            Cell[i][j].A = area(i, j);
        } /* endfor */
    } /* endfor */

    if (BndNorthSpline.np != 0 ) 
       BndNorthSpline.Rotate_Spline(Angle);
    if (BndSouthSpline.np != 0 ) 
       BndSouthSpline.Rotate_Spline(Angle);
    if (BndEastSpline.np != 0 ) 
       BndEastSpline.Rotate_Spline(Angle);
    if (BndWestSpline.np != 0 )
       BndWestSpline.Rotate_Spline(Angle);
 
}


/*!
 * Re-computes the locations of the nodes and cells of  
 * the quadrilateral grid block based on a mirror       
 * reflection about the y=0 axis.  The cells and nodes  
 * are also re-ordered in the i-direction.              
 *                                                      
 */
void Grid2D_Quad_Block_HO::Reflect_Quad_Block(void) {

  int i, j;
  Vector2D *X;
  Spline2D_HO S;

  X = new Vector2D[NNi];

  for ( j = JNl-Nghost ; j <= JNu+Nghost; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      X[NNi-1-i].x = Node[i][j].X.x;
      X[NNi-1-i].y = - Node[i][j].X.y;
    } /* endfor */
    for ( i = INl-Nghost ; i <= INu+Nghost; ++i ) {
      Node[i][j].X = X[i];
    } /* endfor */
  } /* endfor */

  delete []X;
  X = NULL;

  for ( j = JCl-Nghost ; j <= JCu+Nghost ; ++j) {
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
      Cell[i][j].Xc = centroid(i, j);
      Cell[i][j].A = area(i, j);
    } /* endfor */
  } /* endfor */

  if (BndNorthSpline.np != 0 ) 
    BndNorthSpline.Reflect_Spline();
  if (BndSouthSpline.np != 0 ) 
    BndSouthSpline.Reflect_Spline();
  if (BndEastSpline.np != 0 ) 
    BndEastSpline.Reflect_Spline();
  if (BndWestSpline.np != 0 ) {
    BndWestSpline.Reflect_Spline();
    S = BndEastSpline;
    BndEastSpline = BndWestSpline;
    BndWestSpline = S;
    if (S.np != 0) S.deallocate();
  }/* endif */

  /* Require update of the whole mesh */
  Schedule_Interior_Mesh_Update();
  Schedule_Ghost_Cells_Update();

  /* Reset the boundary condition types at the quadrilateral 
     grid block boundaries. */
  
  Set_BCs();

  /* Compute the cells for the quadrilateral mesh block. */

  Update_Cells();
}


/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.                      
 *                                                      
 */
void Grid2D_Quad_Block_HO::Output_Tecplot(const int Block_Number,
					  const int Output_Title,
					  ostream &Out_File) const {
  
  int i, j;
  
  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 << " \\ \n"
	     << "J = " << JNu - JNl + 1 << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 << " \\ \n"
	     << "J = " << JNu - JNl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for (j  = JNl ; j <= JNu ; ++j ) {
    for ( i = INl ; i <= INu ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}


/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.  Includes boundary   
 * nodes.                                               
 */
void Grid2D_Quad_Block_HO::Output_Nodes_Tecplot(const int Block_Number,
						const int Output_Title,
						ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
	     << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << INu - INl + 1 + 2*Nghost << " \\ \n"
	     << "J = " << JNu - JNl + 1 + 2*Nghost << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for (j  = JNl-Nghost ; j <= JNu+Nghost ; ++j ) {
    for ( i = INl-Nghost ; i <= INu+Nghost ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/*!
 * Writes the cells of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with TECPLOT.                      
 */
void Grid2D_Quad_Block_HO::Output_Cells_Tecplot(const int Block_Number,
						const int Output_Title,
						ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "TITLE = \"" << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Cell Locations)"
	     << "\"" << "\n"
	     << "VARIABLES = \"x\" \\ \n"
	     << "\"y\" \n"
	     << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << ICu - ICl + 1 << " \\ \n"
	     << "J = " << JCu - JCl + 1 << " \\ \n"
	     << "F = POINT \n";

  } else {
    Out_File << "ZONE T =  \"Block Number = " << Block_Number
	     << "\" \\ \n"
	     << "I = " << ICu - ICl + 1 << " \\ \n"
	     << "J = " << JCu - JCl + 1 << " \\ \n"
	     << "F = POINT \n";
  } /* endif */

  for (j  = JCl ; j <= JCu ; ++j ) {
    for ( i = ICl ; i <= ICu ; ++i ) {
      Out_File << " " << Cell[i][j].Xc << "\n";
    } /* endfor */
  } /* endfor */
  Out_File << setprecision(6);

}

/*!
 * Writes the nodes of the quadrilateral mesh to the    
 * specified output stream in a format suitable for     
 * plotting the grid with GNUPLOT.                      
 */
void Grid2D_Quad_Block_HO::Output_Gnuplot(const int Block_Number,
					  const int Output_Title,
					  ostream &Out_File) const {

  int i, j;

  Out_File << setprecision(14);
  if (Output_Title) {
    Out_File << "# " << CFFC_Name()
	     << ": 2D Structured Curvilinear Grid Block (Node Locations)"
	     << "\n"
	     << "# x(m), y(m)\n";
  } /* endif */

  for (j  = JNl ; j <= JNu ; ++j ) {
    for ( i = INl ; i <= INu ; ++i ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
    Out_File << "\n";
  } /* endfor */
  for (i  = INl ; i <= INu ; ++i ) {
    for ( j = JNl ; j <= JNu ; ++j ) {
      Out_File << " " << Node[i][j].X << "\n";
    } /* endfor */
    Out_File << "\n";
  } /* endfor */
  Out_File << setprecision(6);

}



/*!
 * Returns a new quadrilateral mesh block with twice    
 * the mesh resolution of the input grid block.         
 *                                                      
 */
void Grid2D_Quad_Block_HO::Double_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original) {

  int i, j, double_resolution_permitted;
  double sp_l, sp_r, sp_m, ds_ratio;
 
  /* Allocate memory for the cells and nodes of the 
     quadrilateral mesh block with twice the resolution. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) { 
    double_resolution_permitted = 0;
  } else {
    double_resolution_permitted = 1;
    allocate(2*(Grid_Original.NCi-2*Grid_Original.Nghost), 
	     2*(Grid_Original.NCj-2*Grid_Original.Nghost),
	     Grid_Original.Nghost);
  } /* endif */

    /* Copy boundary spline info to quadrilateral mesh block 
       with twice the resolution. */

  if (double_resolution_permitted) {

    if (Grid_Original.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original.BndNorthSpline;
    } /* endif */

    if (Grid_Original.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original.BndSouthSpline;
    } /* endif */

    if (Grid_Original.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original.BndEastSpline;
    } /* endif */
  
    if (Grid_Original.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original.BndWestSpline;
    } /* endif */

    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with twice the resolution. */

    SminN = Grid_Original.SminN;
    SmaxN = Grid_Original.SmaxN;
    SminS = Grid_Original.SminS;
    SmaxS = Grid_Original.SmaxS;
    SminE = Grid_Original.SminE;
    SmaxE = Grid_Original.SmaxE;
    SminW = Grid_Original.SminW;
    SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with twice the resolution. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    OrthogonalN = Grid_Original.OrthogonalN;
    OrthogonalS = Grid_Original.OrthogonalS;
    OrthogonalE = Grid_Original.OrthogonalE;
    OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with twice the resolution. */

    for (j  = Grid_Original.JCl ; j <= Grid_Original.JCu ; ++j ) {
      for ( i = Grid_Original.ICl ; i <= Grid_Original.ICu ; ++i ) {
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	  = Grid_Original.nodeSW(i, j).X;
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	  = Grid_Original.xfaceS(i, j);
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	  = Grid_Original.xfaceW(i, j);
	Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	  = Grid_Original.Cell[i][j].Xc;
	if (j == Grid_Original.JCu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.nodeNW(i, j).X;
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.xfaceN(i, j);
	} /* endif */
	if (i == Grid_Original.ICu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
	    = Grid_Original.nodeSE(i, j).X;
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
	    = Grid_Original.xfaceE(i, j);
	} /* endif */
	if (i == Grid_Original.ICu && j == Grid_Original.JCu) {
	  Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
	    [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
	    = Grid_Original.nodeNE(i, j).X;
	} /* endif */
      } /* endfor */
    } /* endfor */

    if (BndWestSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INl][j-1].X, 
		    BndWestSpline);
	sp_r = getS(Node[INl][j+1].X, 
		    BndWestSpline);
	ds_ratio = abs(Node[INl+1][j].X-
		       Node[INl+1][j-1].X)/
	  abs(Node[INl+1][j+1].X-
	      Node[INl+1][j-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INl][j].X = 
	  Spline(sp_m, BndWestSpline);
      } /* endfor */
    } /* endif */

    if (BndEastSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INu][j-1].X, 
		    BndEastSpline);
	sp_r = getS(Node[INu][j+1].X, 
		    BndEastSpline);
	ds_ratio = abs(Node[INu-1][j].X-
		       Node[INu-1][j-1].X)/
	  abs(Node[INu-1][j+1].X-
	      Node[INu-1][j-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INu][j].X = 
	  Spline(sp_m, BndEastSpline);
      } /* endfor */
    } /* endif */

    if (BndSouthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNl].X, 
		    BndSouthSpline);
	sp_r = getS(Node[i+1][JNl].X, 
		    BndSouthSpline);
	ds_ratio = abs(Node[i][JNl+1].X-
		       Node[i-1][JNl+1].X)/
	  abs(Node[i+1][JNl+1].X-
	      Node[i-1][JNl+1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNl].X =  
	  Spline(sp_m, BndSouthSpline);
      } /* endfor */
    } /* endif */

    if (BndNorthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNu].X, 
		    BndNorthSpline);
	sp_r = getS(Node[i+1][JNu].X, 
		    BndNorthSpline);
	ds_ratio = abs(Node[i][JNu-1].X-
		       Node[i-1][JNu-1].X)/
	  abs(Node[i+1][JNu-1].X-
	      Node[i-1][JNu-1].X);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNu].X =  
	  Spline(sp_m, BndNorthSpline);
      } /* endfor */
    } /* endif */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for quadrilateral mesh block 
       with twice the resolution. */

    Set_BCs();

    /* Compute the exterior nodes for quadrilateral mesh block 
       with twice the resolution. */

    Update_Exterior_Nodes();

    /* Compute the cells for quadrilateral mesh block 
       with twice the resolution. */

    Update_Cells();

  } /* endif */

}

/*!
 * Returns a new quadrilateral mesh block with half the 
 * mesh resolution of the input grid block.             
 *                                                      
 */
void Grid2D_Quad_Block_HO::Half_Mesh_Resolution(const Grid2D_Quad_Block_HO &Grid_Original) {

  int i, j, half_resolution_permitted;
 
  /* Allocate memory for the cells and nodes of the 
     quadrilateral mesh block with half the resolution. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) {
    half_resolution_permitted = 0;
  } else {
    half_resolution_permitted = 1;
    allocate((Grid_Original.NCi-2*Grid_Original.Nghost)/2, 
	     (Grid_Original.NCj-2*Grid_Original.Nghost)/2,
	     Grid_Original.Nghost);
  } /* endif */

    /* Copy boundary spline info to quadrilateral mesh block 
       with half the resolution. */

  if (half_resolution_permitted) {

    if (Grid_Original.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original.BndNorthSpline;
    } /* endif */

    if (Grid_Original.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original.BndSouthSpline;
    } /* endif */

    if (Grid_Original.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original.BndEastSpline;
    } /* endif */
  
    if (Grid_Original.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original.BndWestSpline;
    } /* endif */

    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with half the resolution. */

    SminN = Grid_Original.SminN;
    SmaxN = Grid_Original.SmaxN;
    SminS = Grid_Original.SminS;
    SmaxS = Grid_Original.SmaxS;
    SminE = Grid_Original.SminE;
    SmaxE = Grid_Original.SmaxE;
    SminW = Grid_Original.SminW;
    SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with half the resolution. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    OrthogonalN = Grid_Original.OrthogonalN;
    OrthogonalS = Grid_Original.OrthogonalS;
    OrthogonalE = Grid_Original.OrthogonalE;
    OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with half the resolution. */

    for (j  = JNl ; j <= JNu ; ++j ) {
      for ( i = INl ; i <= INu ; ++i ) {
	Node[i][j].X = 
	  Grid_Original.Node[2*(i-INl)+INl]
	  [2*(j-JNl)+JNl].X;
      } /* endfor */
    } /* endfor */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for quadrilateral mesh block 
       with half the resolution. */

    Set_BCs();

    /* Compute the exterior nodes for quadrilateral mesh block 
       with half the resolution. */

    Update_Exterior_Nodes();

    /* Compute the cells for quadrilateral mesh block 
       with half the resolution. */

    Update_Cells();

  } /* endif */

}

/*!
 * Returns a new quadrilateral mesh block containing    
 * one of the specified sectors of the original grid    
 * block with twice the mesh resolution.                
 *                                                      
 */
void Grid2D_Quad_Block_HO::Refine_Mesh(const Grid2D_Quad_Block_HO &Grid_Original,
				       const int Sector) {

  int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;

  /* Allocate memory for the cells and nodes for the 
     refined quadrilateral mesh block. */

  if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
       (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
       (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
       (Grid_Original.Node == NULL) ) {
    mesh_refinement_permitted = 0;
  } else {
    mesh_refinement_permitted = 1;
    allocate(Grid_Original.NCi-2*Grid_Original.Nghost, 
	     Grid_Original.NCj-2*Grid_Original.Nghost,
	     Grid_Original.Nghost);
  } /* endif */

    /* Copy boundary spline info for the refined
       quadrilateral mesh block. */

  if (mesh_refinement_permitted) {

    if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NW ||
	 Sector == GRID2D_QUAD_BLOCK_SECTOR_NE) &&
	Grid_Original.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original.BndNorthSpline;
    } /* endif */

    if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_SE ||
	 Sector == GRID2D_QUAD_BLOCK_SECTOR_SW) &&
	Grid_Original.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original.BndSouthSpline;
    } /* endif */

    if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NE ||
	 Sector == GRID2D_QUAD_BLOCK_SECTOR_SE) &&
	Grid_Original.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original.BndEastSpline;
    } /* endif */
  
    if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NW ||
	 Sector == GRID2D_QUAD_BLOCK_SECTOR_SW) &&
	Grid_Original.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original.BndWestSpline;
    } /* endif */

    /* Assign boundary spline pathlength info for the refined
       quadrilateral mesh block. */

    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      if (Grid_Original.BndNorthSpline.np != 0) {
	SminN = Grid_Original.SminN;
	SmaxN = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNu].X, 
		     Grid_Original.BndNorthSpline);
      } else {
	SminN = ZERO;
	SmaxN = ZERO;
      } /* endif */
      SminS = ZERO;
      SmaxS = ZERO;
      SminE = ZERO;
      SmaxE = ZERO;
      if (Grid_Original.BndWestSpline.np != 0) {
	SminW = getS(Grid_Original.Node[Grid_Original.INl]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndWestSpline);
	SmaxW = Grid_Original.SmaxW;
      } else {
	SminW = ZERO;
	SmaxW = ZERO;
      } /* endif */
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      if (Grid_Original.BndNorthSpline.np != 0) {
	SminN = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNu].X, 
		     Grid_Original.BndNorthSpline);
	SmaxN = Grid_Original.SmaxN;
      } else {
	SminN = ZERO;
	SmaxN = ZERO;
      } /* endif */
      SminS = ZERO;
      SmaxS = ZERO;
      if (Grid_Original.BndEastSpline.np != 0) {
	SminE = getS(Grid_Original.Node[Grid_Original.INu]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndEastSpline);
	SmaxE = Grid_Original.SmaxE;
      } else {
	SminE = ZERO;
	SmaxE = ZERO;
      } /* endif */
      SminW = ZERO;
      SmaxW = ZERO;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      SminN = ZERO;
      SmaxN = ZERO;
      if (Grid_Original.BndSouthSpline.np != 0) {
	SminS = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNl].X, 
		     Grid_Original.BndSouthSpline);
	SmaxS = Grid_Original.SmaxS;
      } else {
	SminS = ZERO;
	SmaxS = ZERO;
      } /* endif */
      if (Grid_Original.BndEastSpline.np != 0) {
	SminE = Grid_Original.SminE;
	SmaxE = getS(Grid_Original.Node[Grid_Original.INu]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndEastSpline);
	SmaxE = Grid_Original.SmaxE;
      } else {
	SminE = ZERO;
	SmaxE = ZERO;
      } /* endif */
      SminW = ZERO;
      SmaxW = ZERO;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      SminN = ZERO;
      SmaxN = ZERO;
      if (Grid_Original.BndSouthSpline.np != 0) {
	SminS = Grid_Original.SminS;
	SmaxS = getS(Grid_Original.Node[Grid_Original.INl+
					(Grid_Original.INu-
					 Grid_Original.INl)/2]
		     [Grid_Original.JNl].X, 
		     Grid_Original.BndSouthSpline);
      } else {
	SminS = ZERO;
	SmaxS = ZERO;
      } /* endif */
      SminE = ZERO;
      SmaxE = ZERO;
      if (Grid_Original.BndWestSpline.np != 0) {
	SminW = Grid_Original.SminW;
	SmaxW = getS(Grid_Original.Node[Grid_Original.INl]
		     [Grid_Original.JNl+
		      (Grid_Original.JNu-
		       Grid_Original.JNl)/2].X, 
		     Grid_Original.BndWestSpline);
      } else {
	SminW = ZERO;
	SmaxW = ZERO;
      } /* endif */
      break;
    } /* endswitch */

    /* Copy node stretching info to refined quadrilateral 
       mesh block. */

    StretchI = Grid_Original.StretchI;
    BetaI = Grid_Original.BetaI;
    TauI = Grid_Original.TauI;
    StretchJ = Grid_Original.StretchJ;
    BetaJ = Grid_Original.BetaJ;
    TauJ = Grid_Original.TauJ;
    if (Grid_Original.BndNorthSpline.np != 0) {
      OrthogonalN = Grid_Original.OrthogonalN;
    } else {
      OrthogonalN = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndSouthSpline.np != 0) {
      OrthogonalS = Grid_Original.OrthogonalS;
    } else {
      OrthogonalS = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndEastSpline.np != 0) {
      OrthogonalE = Grid_Original.OrthogonalE;
    } else {
      OrthogonalE = ORTHOGONAL;
    } /* endif */
    if (Grid_Original.BndWestSpline.np != 0) {
      OrthogonalW = Grid_Original.OrthogonalW;
    } else {
      OrthogonalW = ORTHOGONAL;
    } /* endif */

    // Force orthogonality at all refined mesh boundaries.
    OrthogonalN = ORTHOGONAL;
    OrthogonalS = ORTHOGONAL;
    OrthogonalE = ORTHOGONAL;
    OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for refined 
       quadrilateral mesh block. */

    switch(Sector) {
    case GRID2D_QUAD_BLOCK_SECTOR_NW :
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_NE :
      i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
      i_max = Grid_Original.ICu;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SE :
      i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
      i_max = Grid_Original.ICu;
      j_min = Grid_Original.JCl; 
      j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
      break;
    case GRID2D_QUAD_BLOCK_SECTOR_SW :
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl; 
      j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
      break;
    default:
      i_min = Grid_Original.ICl;
      i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
      j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
      j_max = Grid_Original.JCu;
      break;
    } /* endswitch */

    for ( j  = j_min; j <= j_max ; ++j ) {
      for ( i = i_min ; i <= i_max ; ++i ) {
	Node[2*(i-i_min)+INl  ]
	  [2*(j-j_min)+JNl  ].X 
	  = Grid_Original.nodeSW(i, j).X;
	Node[2*(i-i_min)+INl+1]
	  [2*(j-j_min)+JNl  ].X 
	  = Grid_Original.xfaceS(i, j);
	Node[2*(i-i_min)+INl  ]
	  [2*(j-j_min)+JNl+1].X 
	  = Grid_Original.xfaceW(i, j);
	//     	       Node[2*(i-i_min)+INl+1]
	//                              [2*(j-j_min)+JNl+1].X 
	//                   = Grid_Original.Cell[i][j].Xc;
	Node[2*(i-i_min)+INl+1]
	  [2*(j-j_min)+JNl+1].X 
	  = (Grid_Original.nodeSW(i,j).X +
	     Grid_Original.nodeSE(i,j).X +
	     Grid_Original.nodeNW(i,j).X +
	     Grid_Original.nodeNE(i,j).X)/FOUR;
	if (j == j_max) {
	  Node[2*(i-i_min)+INl  ]
	    [2*(j-j_min)+JNl+2].X 
	    = Grid_Original.nodeNW(i, j).X;
	  Node[2*(i-i_min)+INl+1]
	    [2*(j-j_min)+JNl+2].X 
	    = Grid_Original.xfaceN(i, j);
	} /* endif */
	if (i == i_max) {
	  Node[2*(i-i_min)+INl+2]
	    [2*(j-j_min)+JNl  ].X 
	    = Grid_Original.nodeSE(i, j).X;
	  Node[2*(i-i_min)+INl+2]
	    [2*(j-j_min)+JNl+1].X 
	    = Grid_Original.xfaceE(i, j);
	} /* endif */
	if (i == i_max && j == j_max) {
	  Node[2*(i-i_min)+INl+2]
	    [2*(j-j_min)+JNl+2].X 
	    = Grid_Original.nodeNE(i, j).X;
	} /* endif */
      } /* endfor */
    } /* endfor */

    if (BndWestSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INl][j-1].X, 
		    BndWestSpline);
	sp_r = getS(Node[INl][j+1].X, 
		    BndWestSpline);
	dl = abs(Node[INl][j  ].X - 
		 Node[INl][j-1].X);
	dr = abs(Node[INl][j+1].X - 
		 Node[INl][j  ].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INl][j].X = 
	  Spline(sp_m, BndWestSpline);
      } /* endfor */
    } /* endif */

    if (BndEastSpline.np != 0) {
      for (j  = JNl+1 ; j < JNu ; j += 2 ) {
	sp_l = getS(Node[INu][j-1].X, 
		    BndEastSpline);
	sp_r = getS(Node[INu][j+1].X, 
		    BndEastSpline);
	dl = abs(Node[INu][j  ].X - 
		 Node[INu][j-1].X);
	dr = abs(Node[INu][j+1].X - 
		 Node[INu][j  ].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[INu][j].X = 
	  Spline(sp_m, BndEastSpline);
      } /* endfor */
    } /* endif */

    if (BndSouthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNl].X, 
		    BndSouthSpline);
	sp_r = getS(Node[i+1][JNl].X, 
		    BndSouthSpline);
	dl = abs(Node[i  ][JNl].X - 
		 Node[i-1][JNl].X);
	dr = abs(Node[i+1][JNl].X - 
		 Node[i  ][JNl].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNl].X =  
	  Spline(sp_m, BndSouthSpline);
      } /* endfor */
    } /* endif */

    if (BndNorthSpline.np != 0) {
      for ( i = INl+1 ; i < INu ; i += 2 ) {
	sp_l = getS(Node[i-1][JNu].X, 
		    BndNorthSpline);
	sp_r = getS(Node[i+1][JNu].X, 
		    BndNorthSpline);
	dl = abs(Node[i  ][JNu].X - 
		 Node[i-1][JNu].X);
	dr = abs(Node[i+1][JNu].X - 
		 Node[i  ][JNu].X);
	ds_ratio = dl/(dl+dr);
	sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	Node[i][JNu].X =  
	  Spline(sp_m, BndNorthSpline);
      } /* endfor */
      // 	 for ( i = ICl; i <= ICu; i++) {
      // 	   if (area(i,JCu) <= ZERO) {
      // 	     Node[i][JNu-1].X = Node[i][JNu-3].X + (2.0/3.0)*(Node[i][JNu].X-
      // 												      Node[i][JNu-3].X);
      // 	     Node[i][JNu-2].X = Node[i][JNu-3].X + (1.0/3.0)*(Node[i][JNu].X-
      // 												      Node[i][JNu-3].X);
      // 	     Node[i+1][JNu-1].X = Node[i+1][JNu-3].X + (2.0/3.0)*(Node[i+1][JNu].X-
      // 												      Node[i+1][JNu-3].X);
      // 	     Node[i+1][JNu-2].X = Node[i+1][JNu-3].X + (1.0/3.0)*(Node[i+1][JNu].X-
      // 												      Node[i+1][JNu-3].X);
      // 	   }
      // 	 }
    } /* endif */

    /* Require update of the interior cells geometric properties. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for refined 
       quadrilateral mesh block. */

    Set_BCs();

    /* Compute the exterior nodes for refined 
       quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for refined 
       quadrilateral mesh block. */

    Update_Cells();

  } /* endif */

}

/********************************************************
 * Routine: Coarsen_Mesh                                *
 *                                                      *
 * Returns a new quadrilateral mesh block resulting     *
 * from the coarsening of four original grid blocks     *
 * with half the original mesh resolution.              *
 *                                                      *
 ********************************************************/
void Grid2D_Quad_Block_HO::Coarsen_Mesh(const Grid2D_Quad_Block_HO &Grid_Original_SW,
					const Grid2D_Quad_Block_HO &Grid_Original_SE,
					const Grid2D_Quad_Block_HO &Grid_Original_NW,
					const Grid2D_Quad_Block_HO &Grid_Original_NE) {

  int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;
 
  /* Allocate memory for the cells and nodes for the 
     coarsened quadrilateral mesh block. */

  if ( (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost-
	2*((Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost)/2) != 0) || 
       (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost-
	2*((Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost)/2) != 0) ||
       (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
       (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
       (Grid_Original_SE.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_SE.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_NW.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_NW.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_NE.NCi != Grid_Original_SW.NCi) ||
       (Grid_Original_NE.NCj != Grid_Original_SW.NCj) ||
       (Grid_Original_SW.Node == NULL) ||
       (Grid_Original_SE.Node == NULL) ||
       (Grid_Original_NW.Node == NULL) ||
       (Grid_Original_NE.Node == NULL) ) {
    mesh_coarsening_permitted = 0;
  } else {
    mesh_coarsening_permitted = 1;
    allocate(Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost, 
	     Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost,
	     Grid_Original_SW.Nghost);
  } /* endif */

    /* Copy boundary spline info for the coarsened
       quadrilateral mesh block. */

  if (mesh_coarsening_permitted) {

    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      BndNorthSpline = Grid_Original_NW.BndNorthSpline;
    } /* endif */

    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      BndSouthSpline = Grid_Original_SW.BndSouthSpline;
    } /* endif */

    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      BndEastSpline = Grid_Original_SE.BndEastSpline;
    } /* endif */
  
    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      BndWestSpline = Grid_Original_SW.BndWestSpline;
    } /* endif */

    /* Assign boundary spline pathlength info for the coarsened
       quadrilateral mesh block. */

    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      SminN = Grid_Original_NW.SminN;
      SmaxN = Grid_Original_NE.SmaxN;
    } else {
      SminN = ZERO;
      SmaxN = ZERO;
    } /* endif */

    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      SminS = Grid_Original_SW.SminS;
      SmaxS = Grid_Original_SE.SmaxS;
    } else {
      SminS = ZERO;
      SmaxS = ZERO;
    } /* endif */

    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      SminE = Grid_Original_SE.SmaxE;
      SmaxE = Grid_Original_NE.SmaxE;
    } else {
      SminE = ZERO;
      SmaxE = ZERO;
    } /* endif */

    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      SminW = Grid_Original_SW.SminW;
      SmaxW = Grid_Original_NW.SmaxW;
    } else {
      SminW = ZERO;
      SmaxW = ZERO;
    } /* endif */

    /* Copy node stretching info to coarsened quadrilateral 
       mesh block. */

    StretchI = Grid_Original_SW.StretchI;
    BetaI = Grid_Original_SW.BetaI;
    TauI = Grid_Original_SW.TauI;
    StretchJ = Grid_Original_SW.StretchJ;
    BetaJ = Grid_Original_SW.BetaJ;
    TauJ = Grid_Original_SW.TauJ;
    if (Grid_Original_NW.BndNorthSpline.np != 0 &&
	Grid_Original_NE.BndNorthSpline.np != 0) {
      OrthogonalN = Grid_Original_NW.OrthogonalN;
    } else {
      OrthogonalN = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SW.BndSouthSpline.np != 0 &&
	Grid_Original_SE.BndSouthSpline.np != 0) {
      OrthogonalS = Grid_Original_SW.OrthogonalS;
    } else {
      OrthogonalS = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SE.BndEastSpline.np != 0 &&
	Grid_Original_NE.BndEastSpline.np != 0) {
      OrthogonalE = Grid_Original_SE.OrthogonalE;
    } else {
      OrthogonalE = ORTHOGONAL;
    } /* endif */
    if (Grid_Original_SW.BndWestSpline.np != 0 &&
	Grid_Original_NW.BndWestSpline.np != 0) {
      OrthogonalW = Grid_Original_SW.OrthogonalW;
    } else {
      OrthogonalW = ORTHOGONAL;
    } /* endif */

    // Force orthogonality at all coarsened mesh boundaries.
    OrthogonalN = ORTHOGONAL;
    OrthogonalS = ORTHOGONAL;
    OrthogonalE = ORTHOGONAL;
    OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for the coarsened
       quadrilateral mesh block. */

    for ( j = Grid_Original_SW.JNl; j <= Grid_Original_SW.JNu ; j += 2 ) {
      for ( i = Grid_Original_SW.INl ; i <= Grid_Original_SW.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_SW.INl)/2+
	  INl;
	j_coarse = (j-Grid_Original_SW.JNl)/2+
	  JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_SW.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_SE.JNl; j <= Grid_Original_SE.JNu ; j += 2 ) {
      for ( i = Grid_Original_SE.INl ; i <= Grid_Original_SE.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_SE.INl)/2+
	  (INu-INl)/2+INl;
	j_coarse = (j-Grid_Original_SE.JNl)/2+
	  JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_SE.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_NW.JNl; j <= Grid_Original_NW.JNu ; j += 2 ) {
      for ( i = Grid_Original_NW.INl ; i <= Grid_Original_NW.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_NW.INl)/2+
	  INl;
	j_coarse = (j-Grid_Original_NW.JNl)/2+
	  (JNu-JNl)/2+JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_NW.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    for ( j = Grid_Original_NE.JNl; j <= Grid_Original_NE.JNu ; j += 2 ) {
      for ( i = Grid_Original_NE.INl ; i <= Grid_Original_NE.INu ; i += 2 ) {
	i_coarse = (i-Grid_Original_NE.INl)/2+
	  (INu-INl)/2+INl;
	j_coarse = (j-Grid_Original_NE.JNl)/2+
	  (JNu-JNl)/2+JNl;
	Node[i_coarse][j_coarse].X = Grid_Original_NE.Node[i][j].X;
      } /* endfor */
    } /* endfor */

    /* Require update of the interior cells geometric properties
       for newly coarsed quadrilateral mesh block. */
    Schedule_Interior_Mesh_Update();

    /* Set the boundary condition types for newly coarsened
       quadrilateral mesh block. */

    Set_BCs();

    /* Compute the exterior nodes for newly coarsened 
       quadrilateral mesh block. */

    Update_Exterior_Nodes();

    /* Compute the cells for newly coarsened
       quadrilateral mesh block. */

    Update_Cells();

  } /* endif */

}


/*!
 * Adjusts the locations of the boundary nodes of a     
 * quadrilateral grid block so that the new node        
 * locations will match with cell volumes of adjacent   
 * quadrilateral mesh blocks that have lower levels of  
 * mesh refinement (i.e., are coarser mesh blocks).     
 */
void Grid2D_Quad_Block_HO::Fix_Refined_Mesh_Boundaries(const int Fix_North_Boundary,
						       const int Fix_South_Boundary,
						       const int Fix_East_Boundary,
						       const int Fix_West_Boundary) {

  int i, j;
  double ds_ratio, dl, dr;
 
  /* Adjust the node locations of the north boundary
     of the quadrilateral mesh block. */

  if (Fix_North_Boundary) {
    for ( i = INl+1 ; i <= INu-1 ; i+=2 ) {
      dl = abs(Node[i  ][JNu].X - 
	       Node[i-1][JNu].X);
      dr = abs(Node[i+1][JNu].X - 
	       Node[i  ][JNu].X);
      ds_ratio = dl/(dl+dr);
      Node[i][JNu].X = 
	Node[i-1][JNu].X +
	ds_ratio*(Node[i+1][JNu].X-
		  Node[i-1][JNu].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the south boundary
       of the quadrilateral mesh block. */

  if (Fix_South_Boundary) {
    for ( i = INl+1 ; i <= INu-1 ; i+=2 ) {
      dl = abs(Node[i  ][JNl].X - 
	       Node[i-1][JNl].X);
      dr = abs(Node[i+1][JNl].X - 
	       Node[i  ][JNl].X);
      ds_ratio = dl/(dl+dr);
      Node[i][JNl].X = 
	Node[i-1][JNl].X +
	ds_ratio*(Node[i+1][JNl].X-
		  Node[i-1][JNl].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the east boundary
       of the quadrilateral mesh block. */

  if (Fix_East_Boundary) {
    for ( j  = JNl+1; j <= JNu-1; j+=2 ) {
      dl = abs(Node[INu][j  ].X - 
	       Node[INu][j-1].X);
      dr = abs(Node[INu][j+1].X - 
	       Node[INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      Node[INu][j].X = 
	Node[INu][j-1].X +
	ds_ratio*(Node[INu][j+1].X-
		  Node[INu][j-1].X);
    } /* endfor */
  } /* endif */

    /* Adjust the node locations of the west boundary
       of the quadrilateral mesh block. */

  if (Fix_West_Boundary) {
    for ( j  = JNl+1; j <= JNu-1; j+=2 ) {
      dl = abs(Node[INl][j  ].X - 
	       Node[INl][j-1].X);
      dr = abs(Node[INl][j+1].X - 
	       Node[INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      Node[INl][j].X = 
	Node[INl][j-1].X +
	ds_ratio*(Node[INl][j+1].X-
		  Node[INl][j-1].X);
    } /* endfor */
  }/* endif */
  
  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Reset the boundary condition types at the grid block
     boundaries. */
 
  Set_BCs();

  /* Recompute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Recompute the cells for the quadrilateral mesh block. */

  Update_Cells();

}


/*!
 * Returns the adjusted the locations of the boundary   
 * nodes of a quadrilateral grid block to their         
 * original unmodified positions.                       
 *                                                      
 */
void Grid2D_Quad_Block_HO::Unfix_Refined_Mesh_Boundaries(void) {

  int i, j;
  double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
 
  /* Return the nodes at the north boundary
     to their original positions. */

  if (BndNorthSpline.np != 0) {
    for ( i = INl+1 ; i < INu ; i += 2 ) {
      sp_l = getS(Node[i-1][JNu].X, 
		  BndNorthSpline);
      sp_r = getS(Node[i+1][JNu].X, 
		  BndNorthSpline);
      dl = abs(Node[i  ][JNu].X - 
	       Node[i-1][JNu].X);
      dr = abs(Node[i+1][JNu].X - 
	       Node[i  ][JNu].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[i][JNu].X = Spline(sp_m, BndNorthSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the south boundary
       to their original positions. */

  if (BndSouthSpline.np != 0) {
    for ( i = INl+1 ; i < INu ; i += 2 ) {
      sp_l = getS(Node[i-1][JNl].X, 
		  BndSouthSpline);
      sp_r = getS(Node[i+1][JNl].X, 
		  BndSouthSpline);
      dl = abs(Node[i  ][JNl].X - 
	       Node[i-1][JNl].X);
      dr = abs(Node[i+1][JNl].X - 
	       Node[i  ][JNl].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[i][JNl].X = Spline(sp_m, BndSouthSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the east boundary
       to their original positions. */

  if (BndEastSpline.np != 0) {
    for (j  = JNl+1 ; j < JNu ; j += 2 ) {
      sp_l = getS(Node[INu][j-1].X, 
		  BndEastSpline);
      sp_r = getS(Node[INu][j+1].X, 
		  BndEastSpline);
      dl = abs(Node[INu][j  ].X - 
	       Node[INu][j-1].X);
      dr = abs(Node[INu][j+1].X - 
	       Node[INu][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[INu][j].X = Spline(sp_m, BndEastSpline);
    } /* endfor */
  } /* endif */

    /* Return the nodes at the west boundary
       to their original positions. */

  if (BndWestSpline.np != 0) {
    for (j  = JNl+1 ; j < JNu ; j += 2 ) {
      sp_l = getS(Node[INl][j-1].X, 
		  BndWestSpline);
      sp_r = getS(Node[INl][j+1].X, 
		  BndWestSpline);
      dl = abs(Node[INl][j  ].X - 
	       Node[INl][j-1].X);
      dr = abs(Node[INl][j+1].X - 
	       Node[INl][j  ].X);
      ds_ratio = dl/(dl+dr);
      sp_m = sp_l + ds_ratio*(sp_r-sp_l);
      Node[INl][j].X = Spline(sp_m, BndWestSpline);
    } /* endfor */
  }/* endif */

  /* Require update of the interior cells geometric properties. */
  Schedule_Interior_Mesh_Update();

  /* Reset the boundary condition types at the grid block
     boundaries. */
 
  Set_BCs();

  /* Recompute the exterior nodes for the quadrilateral mesh block. */

  Update_Exterior_Nodes();

  /* Recompute the cells for the quadrilateral mesh block. */

  Update_Cells();

}

/*!
 * Output operator.
 */
ostream &operator << (ostream &out_file, 
		      const Grid2D_Quad_Block_HO &G) {
  int i, j;
  out_file << G.NNi << " " << G.INl << " " << G.INu << "\n";
  out_file << G.NNj << " " << G.JNl << " " << G.JNu << "\n";
  out_file << G.Nghost << "\n";
  if (G.NNi == 0 || G.NNj == 0) return(out_file);
  out_file << G.NCi << " " << G.ICl << " " << G.ICu << "\n";
  out_file << G.NCj << " " << G.JCl << " " << G.JCu << "\n";

  // Output node data
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
    for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
      out_file << G.Node[i][j].X << "\n";
    } /* endfor */
  } /* endfor */

  for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
    out_file << G.BCtypeN[i] << " " << G.BCtypeS[i] << "\n";
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
    out_file << G.BCtypeE[j] << " " << G.BCtypeW[j] << "\n";
  } /* endfor */

  // Output North boundary spline information
  if (G.BndNorthSpline.np != 0 ) {
    out_file << G.BndNorthSpline;
  } else {
    out_file << G.BndNorthSpline.np << "\n";
  } /* endif */

  // Output South boundary spline information
  if (G.BndSouthSpline.np != 0 ) {
    out_file << G.BndSouthSpline;
  } else {
    out_file << G.BndSouthSpline.np << "\n";
  } /* endif */

  // Output East boundary spline information
  if (G.BndEastSpline.np != 0 ) {
    out_file << G.BndEastSpline;
  } else {
    out_file << G.BndEastSpline.np << "\n";
  } /* endif */

  // Output West boundary spline information
  if (G.BndWestSpline.np != 0 ) {
    out_file << G.BndWestSpline;
  } else {
    out_file << G.BndWestSpline.np << "\n";
  } /* endif */

  out_file.setf(ios::scientific);
  out_file << G.SminN << " " << G.SmaxN << " " << G.SminS << " " << G.SmaxS << "\n"; 
  out_file << G.SminE << " " << G.SmaxE << " " << G.SminW << " " << G.SmaxW << "\n";
  out_file << G.StretchI << " " << G.BetaI << " " << G.TauI << "\n";
  out_file << G.StretchJ << " " << G.BetaJ << " " << G.TauJ << "\n";
  out_file.unsetf(ios::scientific);
  out_file << G.OrthogonalN << " "  << G.OrthogonalS << " "  
           << G.OrthogonalE << " "  << G.OrthogonalW << "\n";
  return (out_file);
}

/*!
 * Input operator.
 */
istream &operator >> (istream &in_file, 
		      Grid2D_Quad_Block_HO &G) {
  int i, j, ni, il, iu, nj, jl, ju, ng;

  // Read mesh parameters
  in_file.setf(ios::skipws);
  // Read indexes for nodes
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file.unsetf(ios::skipws);

  // Provide enough memory for the new mesh
  if (ni == 0 || nj == 0) {
    if (G.Node != NULL) G.deallocate(); return(in_file);
  } /* endif */
  G.allocate(ni-2*ng-1, nj-2*ng-1, ng);

  // Read indexes for cells
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu; in_file >> nj >> jl >> ju;
  in_file.unsetf(ios::skipws);

  // Read the coordinates of the mesh nodes
  for ( j = G.JNl-G.Nghost ; j <= G.JNu+G.Nghost; ++j ) {
    for ( i = G.INl-G.Nghost ; i <= G.INu+G.Nghost; ++i ) {
      in_file >> G.Node[i][j].X;
    } /* endfor */
  } /* endfor */

  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
    for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
      G.Cell[i][j].I = i; G.Cell[i][j].J = j;
      G.Cell[i][j].Xc = G.centroid(i, j);
      G.Cell[i][j].A = G.area(i, j);
    } /* endfor */
  } /* endfor */

  // Read the type of boundary conditions for each boundary
  for ( i = G.ICl-G.Nghost ; i <= G.ICu+G.Nghost ; ++i) {
    in_file.setf(ios::skipws);
    in_file >> G.BCtypeN[i] >> G.BCtypeS[i];
    in_file.unsetf(ios::skipws);
  } /* endfor */
  for ( j = G.JCl-G.Nghost ; j <= G.JCu+G.Nghost ; ++j) {
    in_file.setf(ios::skipws);
    in_file >> G.BCtypeE[j] >> G.BCtypeW[j];
    in_file.unsetf(ios::skipws);
  } /* endfor */

  // Read the North boundary spline
  in_file >> G.BndNorthSpline;

  // Read the South boundary spline
  in_file >> G.BndSouthSpline;

  // Read the East boundary spline
  in_file >> G.BndEastSpline;

  // Read the West boundary spline
  in_file >> G.BndWestSpline;

  in_file.setf(ios::skipws);
  in_file >> G.SminN >> G.SmaxN >> G.SminS >> G.SmaxS; 
  in_file >> G.SminE >> G.SmaxE >> G.SminW >> G.SmaxW;
  in_file >> G.StretchI >> G.BetaI >> G.TauI;
  in_file >> G.StretchJ >> G.BetaJ >> G.TauJ;
  in_file >> G.OrthogonalN >> G.OrthogonalS 
          >> G.OrthogonalE >> G.OrthogonalW;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*
 * Determine the minimum distance between the Node[i][j] 
 * and all the neighbour edges and cell diagonals,
 * for a quadrilateral mesh. 
 * \verbatim
 *
 * For the quadrilateral mesh there are 8 edges and 4 diagonals 
 * that must be checked for determining the min distance.
 *
 *          6      5       
 *       o-----o------o    
 *       |            |    
 *     7 |   (i,j)    | 4  
 *       o     o      o    
 *       |            |    
 *     8 |            | 3  
 *       o-----o------o    
 *          1     2        
 *
 * \endverbatim
 *
 * \param [in] i i-index of the node
 * \param [in] j j-index of the node
 */
double Grid2D_Quad_Block_HO::MinimumNodeEdgeDistance(const int &i, const int &j) const {

  double MinDistance;

  // Edge 1
  MinDistance = DistanceFromPointToLine(Node[i][j].X, Node[i-1][j-1].X, Node[i][j-1].X);
  // Edge 2
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j-1].X, Node[i+1][j-1].X));
  // Edge 3
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j-1].X, Node[i+1][j].X));
  // Edge 4
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j].X, Node[i+1][j+1].X));
  // Edge 5
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j+1].X, Node[i][j+1].X));
  // Edge 6
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j+1].X, Node[i-1][j+1].X));
  // Edge 7
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j+1].X, Node[i-1][j].X));
  // Edge 8
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j].X, Node[i-1][j-1].X));

  // Diagonal 1
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j-1].X, Node[i+1][j].X));
  // Diagonal 2
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i+1][j].X, Node[i][j+1].X));
  // Diagonal 3
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i][j+1].X, Node[i-1][j].X));
  // Diagonal 4
  MinDistance = min(MinDistance, DistanceFromPointToLine(Node[i][j].X, Node[i-1][j].X, Node[i][j-1].X));

  return MinDistance;
}

/*
 * Determines the distance between the point of interest
 * and a line defined by the two points. 
 *
 * \param Point the position vector of the point of interest
 * \param LinePoint1 the position vector of one line end 
 * \param LinePoint2 the position vector of the other line end
 */
double Grid2D_Quad_Block_HO::DistanceFromPointToLine(const Vector2D &Point,
						     const Vector2D &LinePoint1,
						     const Vector2D &LinePoint2) const {

  double DeltaX, DeltaY, X_Dis, Y_Dis;

  // X_Dis, Y_Dis -> the X and Y coordinates of the point which is used to measure the distance

  DeltaX = LinePoint2.x - LinePoint1.x;
  DeltaY = LinePoint2.y - LinePoint1.y;

  // Determine X_Dis
  X_Dis = DeltaX*DeltaX*Point.x + DeltaX*DeltaY*(Point.y - LinePoint1.y) + DeltaY*DeltaY*LinePoint1.x;
  X_Dis /= (DeltaX*DeltaX + DeltaY*DeltaY);

  if ( (fabs(X_Dis - LinePoint1.x) <= 1.0e-14) && (fabs(X_Dis - LinePoint2.x) <= 1.0e-14) ){
    // the line defined by the two points is parallel to Oy
    return fabs(Point.x - LinePoint1.x);
  }

  // Determine Y_Dis
  Y_Dis = LinePoint1.y + DeltaY*(X_Dis - LinePoint1.x)/DeltaX;

  // Return the distance value
  DeltaX = Point.x - X_Dis;
  DeltaY = Point.y - Y_Dis;

  return sqrt( DeltaX*DeltaX + DeltaY*DeltaY );
}

/* \file HO_Grid2DQuad.cc
   \brief Source file initializing/implementing member variables/functions that belong to classes defined in HO_Grid2DQuad.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuad.h"	// Include 2D high-order block grid


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
    SminN(ZERO), SmaxN(ZERO), SminS(ZERO), SmaxS(ZERO), 
    SminE(ZERO), SmaxE(ZERO), SminW(ZERO), SmaxW(ZERO),
    StretchI(0), StretchJ(0), BetaI(ONE), TauI(ONE),
    BetaJ(ONE), TauJ(ONE),
    OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1)
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
   SminN(ZERO), SmaxN(ZERO), SminS(ZERO), SmaxS(ZERO), 
   SminE(ZERO), SmaxE(ZERO), SminW(ZERO), SmaxW(ZERO),
   StretchI(0), StretchJ(0),
   BetaI(ONE), TauI(ONE),
   BetaJ(ONE), TauJ(ONE),
   OrthogonalN(1), OrthogonalS(1), OrthogonalE(1), OrthogonalW(1)
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
	Cell[i][j].I  = G.Cell[i][j].I;
	Cell[i][j].J  = G.Cell[i][j].J;
	Cell[i][j].Xc = G.Cell[i][j].Xc;
	Cell[i][j].A  = G.Cell[i][j].A;
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
    Copy_Spline(BndNorthSpline, G.BndNorthSpline);
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
  } /* endif */

  if (G.BndSouthSpline.np != 0) {
    Copy_Spline(BndSouthSpline, G.BndSouthSpline);
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
  } /* endif */
  
  if (G.BndEastSpline.np != 0) {
    Copy_Spline(BndEastSpline, G.BndEastSpline);
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
  } /* endif */
  
  if (G.BndWestSpline.np != 0) {
    Copy_Spline(BndWestSpline, G.BndWestSpline);
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
}

/*!
 * Allocate memory.
 *
 * \param Ni number of cells in i-direction
 * \param Nj number of cells in j-direction
 * \param Ng number of ghost cells
 */
void Grid2D_Quad_Block_HO::allocate(const int Ni, const int Nj, const int Ng) {
  int i;

  // Check conditions
  assert( Ni > 1 && Nj > 1 && Ng >= 2);

  // Check if the new required memory has dimensions different than the currently allocated ones
  if ( Ni != (NCi-2*Nghost)  ||  Nj != (NCj-2*Nghost) ){

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
    BCtypeN = new int[NCi]; BCtypeS = new int[NCi];
    BCtypeE = new int[NCj]; BCtypeW = new int[NCj];
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
	Cell[i][j].I  = Grid.Cell[i][j].I;
	Cell[i][j].J  = Grid.Cell[i][j].J;
	Cell[i][j].Xc = Grid.Cell[i][j].Xc;
	Cell[i][j].A  = Grid.Cell[i][j].A;
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
    Copy_Spline(BndNorthSpline, Grid.BndNorthSpline);
  } else if (BndNorthSpline.np != 0) {
    BndNorthSpline.deallocate();
  } /* endif */

  if (Grid.BndSouthSpline.np != 0) {
    Copy_Spline(BndSouthSpline, Grid.BndSouthSpline);
  } else if (BndSouthSpline.np != 0) {
    BndSouthSpline.deallocate();
  } /* endif */
  
  if (Grid.BndEastSpline.np != 0) {
    Copy_Spline(BndEastSpline, Grid.BndEastSpline);
  } else if (BndEastSpline.np != 0) {
    BndEastSpline.deallocate();
  } /* endif */
  
  if (Grid.BndWestSpline.np != 0) {
    Copy_Spline(BndWestSpline, Grid.BndWestSpline);
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
 * Output operators.
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
    out_file << G.BndNorthSpline.np   << " " 
	     << G.BndNorthSpline.type << "\n"; 
    out_file << G.BndNorthSpline;
  } else {
    out_file << G.BndNorthSpline.np << "\n";
  } /* endif */

  // Output South boundary spline information
  if (G.BndSouthSpline.np != 0 ) {
    out_file << G.BndSouthSpline.np   << " " 
	     << G.BndSouthSpline.type << "\n"; 
    out_file << G.BndSouthSpline;
  } else {
    out_file << G.BndSouthSpline.np << "\n";
  } /* endif */

  // Output East boundary spline information
  if (G.BndEastSpline.np != 0 ) {
    out_file << G.BndEastSpline.np   << " " 
	     << G.BndEastSpline.type << "\n"; 
    out_file << G.BndEastSpline;
  } else {
    out_file << G.BndEastSpline.np << "\n";
  } /* endif */

  // Output West boundary spline information
  if (G.BndWestSpline.np != 0 ) {
    out_file << G.BndWestSpline.np   << " " 
	     << G.BndWestSpline.type << "\n"; 
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
  if (G.Node == NULL || G.Cell == NULL || G.NNi != ni || G.NNj != nj) {
    G.allocate(ni-2*ng-1, nj-2*ng-1, ng);
  } /* endif */

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
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (G.BndNorthSpline.np != 0) G.BndNorthSpline.deallocate(); 
  if (ni != 0 ) {
    G.BndNorthSpline.allocate(ni);
    in_file.setf(ios::skipws);
    in_file >> i; G.BndNorthSpline.settype(i);
    in_file.unsetf(ios::skipws);
    in_file >> G.BndNorthSpline; G.BndNorthSpline.pathlength();
  } /* endif */

  // Read the South boundary spline
  in_file.setf(ios::skipws);
  in_file >> ni;
  in_file.unsetf(ios::skipws);
  if (G.BndSouthSpline.np != 0) G.BndSouthSpline.deallocate(); 
  if (ni != 0 ) {
    G.BndSouthSpline.allocate(ni);
    in_file.setf(ios::skipws);
    in_file >> i; G.BndSouthSpline.settype(i);
    in_file.unsetf(ios::skipws);
    in_file >> G.BndSouthSpline; G.BndSouthSpline.pathlength();
  } /* endif */

  // Read the East boundary spline
  in_file.setf(ios::skipws);
  in_file >> nj;
  in_file.unsetf(ios::skipws);
  if (G.BndEastSpline.np != 0) G.BndEastSpline.deallocate(); 
  if (nj != 0 ) {
    G.BndEastSpline.allocate(nj);
    in_file.setf(ios::skipws);
    in_file >> j; G.BndEastSpline.settype(j);
    in_file.unsetf(ios::skipws);
    in_file >> G.BndEastSpline; G.BndEastSpline.pathlength();
  } /* endif */

  // Read the West boundary spline
  in_file.setf(ios::skipws);
  in_file >> nj;
  in_file.unsetf(ios::skipws);
  if (G.BndWestSpline.np != 0) G.BndWestSpline.deallocate(); 
  if (nj != 0 ) {
    G.BndWestSpline.allocate(nj);
    in_file.setf(ios::skipws);
    in_file >> j; G.BndWestSpline.settype(j);
    in_file.unsetf(ios::skipws);
    in_file >> G.BndWestSpline; G.BndWestSpline.pathlength();
  } /* endif */


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


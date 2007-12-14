/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DdRdU.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This file contains the function definitions for 
///        Jacobian formulations related to Flame2DState.
///        These Jacobians are used for:
///        -> Semi-Implicit (Source only)
///        -> Newton-Krylov GMRES Preconditioner
///           (1st and 2nd order stencil Inviscid, Viscous, and Source)
///
/////////////////////////////////////////////////////////////////////
#include "Flame2DdRdU.h"

// All functions setup to use "Uo" values as point implicit 
// can be used in a multistage scheme and this insures 
// that all Jacobians are taken as dR/dUo.


/////////////////////////////////////////////////////////////////////
/// Semi-Implicit Block Jacobi formulations
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Routine: SemiImplicitBlockJacobi                     *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi(::DenseMatrix &dSdU,
			     Flame2D_Quad_Block &SolnBlk,
			     const Flame2D_pState &Wo,
			     const int &ii, const int &jj){ 
  
  if (SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID) { 

    // declares
    int n = dSdU.get_n(); 
    static ::DenseMatrix dRdW(n,n); dRdW.zero();
    static ::DenseMatrix dWdU(n,n); dWdU.zero();
    
    // Add Source Jacobians (viscous axisymmetric, turbulence)
    SemiImplicitBlockJacobi_dSdW(dRdW,SolnBlk,Wo,ii,jj);
    
    // Transformation Jacobian 
    Wo.dWdU(dWdU);  // <- need to load dih/dic temporary array (done outside)
    dSdU += dRdW*dWdU;
  }

  // Add Source Jacobians (inviscid axisymmetric, chemistry, gravity)
  SemiImplicitBlockJacobi_dSdU(dSdU,SolnBlk,Wo,ii,jj);                 

}


/********************************************************
 * Routine: SemiImplicitBlockJacobi_dSdW                *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi_dSdW(::DenseMatrix &dSdW,
				  Flame2D_Quad_Block &SolnBlk,
				  const Flame2D_pState &Wo,
				  const int &ii, const int &jj){ 
  
  //
  // Viscous Axisymmetric source term jacobian    
  //
  if(SolnBlk.Axisymmetric && SolnBlk.Flow_Type != FLOWTYPE_INVISCID){     
    
    //Cacluate 2nd derivatives  
    double d_dWdx_dW_C,d_dWdy_dW_C;
    d_dWd_dW_Center(d_dWdx_dW_C,d_dWdy_dW_C,SolnBlk,ii, jj);  
    
    // Compute the jacobian
    Wo.dSa_vdW(dSdW,
	       SolnBlk.dWdx[ii][jj],
	       SolnBlk.dWdy[ii][jj],
	       SolnBlk.Grid.Cell[ii][jj].Xc,
	       SolnBlk.Axisymmetric,
	       d_dWdx_dW_C,d_dWdy_dW_C);
  }
  
}


/********************************************************
 * Routine: SemiImplicitBlockJacobi_dSdU                *
 *                                                      *
 * This routine adds all the appropriate source         *    
 * Jacobians based on flow type, for in semi-implicit   *
 * and implicit calculations.                           *
 *                                                      * 
 ********************************************************/ 
void SemiImplicitBlockJacobi_dSdU(::DenseMatrix &dSdU,
				  Flame2D_Quad_Block &SolnBlk,
				  const Flame2D_pState &Wo,
				  const int &ii, const int &jj){
  
  //Add Jacobian for inviscid axisymmetric source terms
  if (SolnBlk.Axisymmetric) {
    Wo.dSa_idU(dSdU,
	       SolnBlk.Grid.Cell[ii][jj].Xc, 
	       SolnBlk.Axisymmetric);
  }
  
  //Add Jacobian for finite-rate chemistry source terms  
  if (Flame2D_pState::isReacting()) Wo.dSwdU(dSdU);  

  //Add Jacobian for gravitational source terms
  if (SolnBlk.Gravity){
    Wo.dSgdU(dSdU);
  } 
  
}



/********************************************************
 * Routine: d_dWd_dW_Diamond                            *
 *                                                      *
 * This routine calculates the 2nd deriavaites          *
 * associated with diamond path and bilinear            *
 * interpolation.                                       *
 *                                                      *
 ********************************************************/
void d_dWd_dW_Diamond(double &d_dWdx_dW, double &d_dWdy_dW, 
		      Flame2D_Quad_Block &SolnBlk, 
		      const double &LEFT, const double &RIGHT, 
		      const int &Orient_cell, const int &Orient_face,  
		      const int &i, const int &j){

  //  double area[4];
  double AREA;
  Vector2D norm[4];
  double  dWnNWdWc, dWnNEdWc,  dWnSWdWc, dWnSEdWc;
 
  switch(Orient_face){
    /*************** NORTH ****************************/
  case NORTH: 
    dWnNWdWc = LEFT;
    dWnNEdWc = RIGHT;

    //  normal vector of the SE side of a diamond 
    norm[0].x = SolnBlk.Grid.nodeNE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[0].y = -( SolnBlk.Grid.nodeNE(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.Cell[i][j+1].Xc.y - SolnBlk.Grid.nodeNE(i,j).X.y;
    norm[1].y = -(SolnBlk.Grid.Cell[i][j+1].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j+1].Xc.y;
    norm[2].y = -(SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i][j+1].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[3].y = -( SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j+1].Xc)));
    
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnNEdWc)* norm[0].x+dWnNEdWc* norm[1].x+ dWnNWdWc* norm[2].x+ (ONE+dWnNWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnNEdWc)* norm[0].y+dWnNEdWc* norm[1].y+ dWnNWdWc* norm[2].y+ (ONE+dWnNWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == NORTH) {
      d_dWdx_dW = HALF*(dWnNEdWc*norm[0].x + (ONE+dWnNEdWc)* norm[1].x + (ONE+dWnNWdWc)* norm[2].x + dWnNWdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnNEdWc*norm[0].y + (ONE+dWnNEdWc)* norm[1].y + (ONE+dWnNWdWc)* norm[2].y + dWnNWdWc*norm[3].y)/AREA;  
    } else if( Orient_cell == EAST || Orient_cell == NORTH_EAST) {
      d_dWdx_dW = HALF*( dWnNEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*( dWnNEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if( Orient_cell == WEST || Orient_cell == NORTH_WEST) {
      d_dWdx_dW = HALF*( dWnNWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnNWdWc*(norm[2].y+norm[3].y))/AREA;  
    }

    break;
    
    /*************** EAST ****************************/
  case EAST:
    dWnNEdWc = LEFT;
    dWnSEdWc = RIGHT; 

    //  normal vector of the SE side of a diamond 
    norm[0].x =  SolnBlk.Grid.Cell[i+1][j].Xc.y - SolnBlk.Grid.nodeSE(i,j).X.y;
    norm[0].y = -(SolnBlk.Grid.Cell[i+1][j].Xc.x - SolnBlk.Grid.nodeSE(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.nodeNE(i,j).X.y -  SolnBlk.Grid.Cell[i+1][j].Xc.y ;
    norm[1].y = -(SolnBlk.Grid.nodeNE(i,j).X.x -  SolnBlk.Grid.Cell[i+1][j].Xc.x );
    //  normal vector of the NW side of a diamond 
    norm[2].x =   SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeNE(i,j).X.y ;
    norm[2].y = -(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeNE(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x = SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[3].y = -(SolnBlk.Grid.nodeSE(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)^
		       (SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i+1][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeNE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
   
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc* norm[0].x+ dWnNEdWc* norm[1].x+ (ONE+ dWnNEdWc)* norm[2].x+ (ONE+dWnSEdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc* norm[0].y+ dWnNEdWc* norm[1].y+ (ONE+ dWnNEdWc)* norm[2].y+ (ONE+dWnSEdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == EAST) {
      d_dWdx_dW = HALF*( (ONE+dWnSEdWc)* norm[0].x + (ONE+dWnNEdWc)*norm[1].x + dWnNEdWc* norm[2].x + dWnSEdWc*norm[3].x)/AREA;
      d_dWdy_dW = HALF*( (ONE+dWnSEdWc)* norm[0].y + (ONE+dWnNEdWc)*norm[1].y + dWnNEdWc* norm[2].y + dWnSEdWc*norm[3].y)/AREA;  
    } else if( Orient_cell == NORTH || Orient_cell == NORTH_EAST) {
      d_dWdx_dW = HALF*( dWnNEdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*( dWnNEdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if( Orient_cell == SOUTH || Orient_cell == SOUTH_EAST) {
      d_dWdx_dW = HALF*( dWnSEdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSEdWc*(norm[0].y+norm[3].y))/AREA;  
    }
    break;

    /*************** SOUTH ****************************/
  case SOUTH:
    dWnSEdWc = LEFT;
    dWnSWdWc = RIGHT;

    //  normal vector of the SE side of a diamond 
    norm[0].x =  SolnBlk.Grid.nodeSE(i,j).X.y - SolnBlk.Grid.Cell[i][j-1].Xc.y;
    norm[0].y = -(SolnBlk.Grid.nodeSE(i,j).X.x - SolnBlk.Grid.Cell[i][j-1].Xc.x  );
    //  normal vector of the NE side of a diamond 
    norm[1].x = SolnBlk.Grid.Cell[i][j].Xc.y -  SolnBlk.Grid.nodeSE(i,j).X.y;
    norm[1].y = -(SolnBlk.Grid.Cell[i][j].Xc.x -  SolnBlk.Grid.nodeSE(i,j).X.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =   SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y ;
    norm[2].y = -(SolnBlk.Grid.nodeSW(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  SolnBlk.Grid.Cell[i][j-1].Xc.y - SolnBlk.Grid.nodeSW(i,j).X.y;
    norm[3].y = -(SolnBlk.Grid.Cell[i][j-1].Xc.x- SolnBlk.Grid.nodeSW(i,j).X.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j-1].Xc)) +
		  fabs((SolnBlk.Grid.nodeSE(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)));
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*(dWnSEdWc* norm[0].x+ (ONE+dWnSEdWc)* norm[1].x+ (ONE+ dWnSWdWc)* norm[2].x+ (dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSEdWc* norm[0].y+ (ONE+dWnSEdWc)* norm[1].y+ (ONE+ dWnSWdWc)* norm[2].y+ (dWnSWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == SOUTH) {
      d_dWdx_dW = HALF*( (ONE+dWnSEdWc)*norm[0].x + dWnSEdWc*norm[1].x + dWnSWdWc*norm[2].x + (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*( (ONE+dWnSEdWc)*norm[0].y + dWnSEdWc*norm[1].y + dWnSWdWc*norm[2].y + (ONE+dWnSWdWc)* norm[3].y)/AREA;  
    } else if( Orient_cell == EAST || Orient_cell == SOUTH_EAST) {
      d_dWdx_dW = HALF*( dWnSEdWc*(norm[0].x+norm[1].x))/AREA;
      d_dWdy_dW = HALF*( dWnSEdWc*(norm[0].y+norm[1].y))/AREA;  
    } else if( Orient_cell == WEST || Orient_cell == SOUTH_WEST) {
      d_dWdx_dW = HALF*( dWnSWdWc*(norm[2].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSWdWc*(norm[2].y+norm[3].y))/AREA;  
    }
    break;
    
    /*************** WEST ****************************/
  case WEST:
    dWnSWdWc = LEFT;
    dWnNWdWc = RIGHT;
    
    //  normal vector of the SE side of a diamond 
    norm[0].x =   SolnBlk.Grid.Cell[i][j].Xc.y - SolnBlk.Grid.nodeSW(i,j).X.y;
    norm[0].y = -(SolnBlk.Grid.Cell[i][j].Xc.x - SolnBlk.Grid.nodeSW(i,j).X.x);
    //  normal vector of the NE side of a diamond 
    norm[1].x =  SolnBlk.Grid.nodeNW(i,j).X.y - SolnBlk.Grid.Cell[i][j].Xc.y;
    norm[1].y = -(SolnBlk.Grid.nodeNW(i,j).X.x - SolnBlk.Grid.Cell[i][j].Xc.x);
    //  normal vector of the NW side of a diamond 
    norm[2].x =  SolnBlk.Grid.Cell[i-1][j].Xc.y - SolnBlk.Grid.nodeNW(i,j).X.y ;
    norm[2].y = -( SolnBlk.Grid.Cell[i-1][j].Xc.x - SolnBlk.Grid.nodeNW(i,j).X.x);
    //  normal vector of the SW side of a diamond 
    norm[3].x =  SolnBlk.Grid.nodeSW(i,j).X.y - SolnBlk.Grid.Cell[i-1][j].Xc.y;
    norm[3].y = -(SolnBlk.Grid.nodeSW(i,j).X.x - SolnBlk.Grid.Cell[i-1][j].Xc.x);
    
    AREA =  HALF*(fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i][j].Xc)) +
		  fabs((SolnBlk.Grid.nodeNW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)^
		       (SolnBlk.Grid.nodeSW(i,j).X-SolnBlk.Grid.Cell[i-1][j].Xc)));
    
    if(Orient_cell == CENTER){
      d_dWdx_dW = HALF*((ONE+dWnSWdWc)* norm[0].x+ (ONE+dWnNWdWc)* norm[1].x+ dWnNWdWc* norm[2].x+ dWnSWdWc* norm[3].x)/AREA;
      d_dWdy_dW = HALF*((ONE+dWnSWdWc)* norm[0].y+ (ONE+dWnNWdWc)* norm[1].y+ dWnNWdWc* norm[2].y+ dWnSWdWc* norm[3].y)/AREA;  
    }  else if( Orient_cell == WEST) {
      d_dWdx_dW = HALF*(dWnSWdWc*norm[0].x + dWnNWdWc*norm[1].x + (ONE+dWnNWdWc)* norm[2].x+ (ONE+dWnSWdWc)* norm[3].x)/AREA;
      d_dWdy_dW = HALF*(dWnSWdWc*norm[0].y + dWnNWdWc*norm[1].y + (ONE+dWnNWdWc)* norm[2].y+ (ONE+dWnSWdWc)* norm[3].y)/AREA;    
    } else if( Orient_cell == NORTH || Orient_cell == NORTH_WEST) {
      d_dWdx_dW = HALF*( dWnNWdWc*(norm[1].x+norm[2].x))/AREA;
      d_dWdy_dW = HALF*( dWnNWdWc*(norm[1].y+norm[2].y))/AREA;  
    } else if( Orient_cell == SOUTH || Orient_cell == SOUTH_WEST) {
      d_dWdx_dW = HALF*( dWnSWdWc*(norm[0].x+norm[3].x))/AREA;
      d_dWdy_dW = HALF*( dWnSWdWc*(norm[0].y+norm[3].y))/AREA;      
    }
    break;

  }

}

/********************************************************
 * Routine: d_dWd_dW_Center                             *
 *                                                      *
 * This routine calculates the 2nd deriavaites          *
 * associated with diamond path and bilinear            *
 * interpolation.                                       *
 *                                                      *
 ********************************************************/
void d_dWd_dW_Center(double &d_dWdx_dW_C, double &d_dWdy_dW_C, 
		     Flame2D_Quad_Block &SolnBlk, 
		     const int &i, const int &j){

  double area[4], d_dWdx_dW[4], d_dWdy_dW[4];

  // area weighted gradients at cell centers, 4 inside triangles
  area[0] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i][j+1].X )^
    (SolnBlk.Grid.xfaceN(i,j)- SolnBlk.Grid.Cell[i][j].Xc);
  area[1] = HALF*(SolnBlk.Grid.Node[i+1][j+1].X - SolnBlk.Grid.Node[i+1][j].X )^
    (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceE(i,j)); 
  area[2] = HALF*(SolnBlk.Grid.Node[i+1][j].X - SolnBlk.Grid.Node[i][j].X )^
    (SolnBlk.Grid.Cell[i][j].Xc - SolnBlk.Grid.xfaceS(i,j) );  
  area[3] = HALF*(SolnBlk.Grid.Node[i][j+1].X - SolnBlk.Grid.Node[i][j].X )^
    ( SolnBlk.Grid.xfaceW(i, j) - SolnBlk.Grid.Cell[i][j].Xc );
  

  //NORTH
  d_dWd_dW_Diamond(d_dWdx_dW[0] ,d_dWdy_dW[0], SolnBlk,
		   SolnBlk.dWn_dWc(i,j+1, NORTH_WEST), SolnBlk.dWn_dWc(i+1, j+1, NORTH_EAST), 
		   CENTER, NORTH, i, j);  
  //EAST
  d_dWd_dW_Diamond(d_dWdx_dW[1] ,d_dWdy_dW[1], SolnBlk,
		   SolnBlk.dWn_dWc(i+1,j+1, NORTH_EAST), SolnBlk.dWn_dWc(i+1, j, SOUTH_EAST), 
		   CENTER, EAST, i, j);
  //SOUTH
  d_dWd_dW_Diamond(d_dWdx_dW[2] ,d_dWdy_dW[2], SolnBlk,
		   SolnBlk.dWn_dWc(i+1,j, SOUTH_EAST), SolnBlk.dWn_dWc(i, j, SOUTH_WEST), 
		   CENTER, SOUTH, i, j);
  //WEST
  d_dWd_dW_Diamond(d_dWdx_dW[3] ,d_dWdy_dW[3], SolnBlk,
		   SolnBlk.dWn_dWc(i,j, SOUTH_WEST), SolnBlk.dWn_dWc(i, j+1, NORTH_WEST), 
		   CENTER, WEST, i, j);
  
  //2nd derivative's at cell center 
  d_dWdx_dW_C = (d_dWdx_dW[0]*area[0] + d_dWdx_dW[1]*area[1] +
		 d_dWdx_dW[2]*area[2] + d_dWdx_dW[3]*area[3])/SolnBlk.Grid.Cell[i][j].A; 
  
  d_dWdy_dW_C = (d_dWdy_dW[0]*area[0]+d_dWdy_dW[1]*area[1] +
		 d_dWdy_dW[2]*area[2] + d_dWdy_dW[3]*area[3])/SolnBlk.Grid.Cell[i][j].A; 
      

}

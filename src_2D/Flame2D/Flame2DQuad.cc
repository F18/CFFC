/**********************************************************************
 * Flame2DQuad.cc: Subroutines for the 2D chem quadrilateral solution  *
 *                block class.                                        *
 **********************************************************************/

// Include 2D Chem quadrilateral solution block class.

#ifndef _FLAME2D_QUAD_INCLUDED
#include "Flame2DQuad.h"
#endif // _FLAME2D_QUAD_INCLUDED

/**********************************************************************
 * Flame2D_Quad_Block -- Create storage for the total number of        *
 *                      variables.                                    *
 **********************************************************************/
int Flame2D_Quad_Block::residual_variable = 1;
int Flame2D_Quad_Block::Number_of_Residual_Norms = 4;
//int Flame2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;

//SNBCK data object
PlanckMean* Flame2D_Quad_Block::PlanckMean_data=NULL;

/////////////////////////////////////////////////////////////////////
/// Reconstruction Methods
/////////////////////////////////////////////////////////////////////


/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j) of   *
 * the computational mesh for the specified             *
 * quadrilateral solution block.  A Green-Gauss         *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
void Flame2D_Quad_Block::Linear_Reconstruction_GreenGauss(const int i, 
							  const int j,
							  const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi_;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    double l_north, l_south, l_east, l_west;
    static Vector2D n_north, n_south, n_east, n_west;
    static Vector2D dXe, dXw, dXn, dXs;
    double W_face;
    
    const int NUM_VAR_FLAME2D = NumVar();
  
   /* Carry out the limited solution reconstruction in
       the specified cell of the computational mesh. */
    
    // Determine the number of neighbouring cells to
    // be used in the reconstruction procedure.  Away from
    // boundaries this 8 neighbours will be used.
    
    /****************************************************************/
    //FOR VISCOUS -> CHANGED TO USE ALL 8
    if (i == ICl-Nghost || i == ICu+Nghost ||
	j == JCl-Nghost || j == JCu+Nghost) {
      n_pts = 0;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j-1;
      i_index[3] = i-1; j_index[3] = j  ;
      i_index[4] = i+1; j_index[4] = j  ;
      i_index[5] = i-1; j_index[5] = j+1;
      i_index[6] = i  ; j_index[6] = j+1;
      i_index[7] = i+1; j_index[7] = j+1;
    }  
    
    /****************************************************************/
    // Perform reconstruction.    
    if (n_pts > 0) {
        // If 8 neighbours are used, apply Green-Gauss reconstruction
      if (n_pts == 8) {
	const Flame2D_State &W_nw = Wnd[i  ][j+1];//WnNW(i, j);
	const Flame2D_State &W_ne = Wnd[i+1][j+1];//WnNE(i, j);
	const Flame2D_State &W_sw = Wnd[i  ][j  ];//WnSW(i, j);
	const Flame2D_State &W_se = Wnd[i+1][j  ];//WnSE(i, j);
	
	l_north = Grid.lfaceN(i, j);
	l_south = Grid.lfaceS(i, j);
	l_east = Grid.lfaceE(i, j);
	l_west = Grid.lfaceW(i, j);
	
	n_north = Grid.nfaceN(i, j);
	n_south = Grid.nfaceS(i, j);
	n_east = Grid.nfaceE(i, j);
	n_west = Grid.nfaceW(i, j);
	
	for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	  W_face = HALF*(W_nw[k]+W_ne[k])*l_north; 
	  dWdx[i][j][k] = W_face*n_north.x;
	  dWdy[i][j][k] = W_face*n_north.y;
	  
	  W_face = HALF*(W_sw[k]+W_se[k])*l_south; 
	  dWdx[i][j][k] += W_face*n_south.x;
	  dWdy[i][j][k] += W_face*n_south.y;
	  
	  W_face = HALF*(W_ne[k]+W_se[k])*l_east; 
	  dWdx[i][j][k] += W_face*n_east.x;
	  dWdy[i][j][k] += W_face*n_east.y;
	  
	  W_face = HALF*(W_nw[k]+W_sw[k])*l_west; 
	  dWdx[i][j][k] += W_face*n_west.x;
	  dWdy[i][j][k] += W_face*n_west.y;
	  
	  dWdx[i][j][k] /= Grid.Cell[i][j].A;
	  dWdy[i][j][k] /= Grid.Cell[i][j].A;
	}
	
      } /* endif */
	   
	// Calculate slope limiters.    
      if (!Freeze_Limiter) {

	dXe = Grid.xfaceE(i, j) - Grid.Cell[i][j].Xc;
	dXw = Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc;
	dXn = Grid.xfaceN(i, j) - Grid.Cell[i][j].Xc;
	dXs = Grid.xfaceS(i, j) - Grid.Cell[i][j].Xc;

	for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
	  u0Min = W[i][j][n];
	  u0Max = u0Min;
	  for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	    u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	    u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
	  } /* endfor */

	  uQuad[0] = W[i][j][n] + dWdx[i][j][n]*dXe.x + dWdy[i][j][n]*dXe.y ;
	  uQuad[1] = W[i][j][n] + dWdx[i][j][n]*dXw.x + dWdy[i][j][n]*dXw.y ;
	  uQuad[2] = W[i][j][n] + dWdx[i][j][n]*dXn.x + dWdy[i][j][n]*dXn.y ;
	  uQuad[3] = W[i][j][n] + dWdx[i][j][n]*dXs.x + dWdy[i][j][n]*dXs.y ;
	    
	  switch(Limiter) {
	  case LIMITER_ONE :
	    phi_ = ONE;
	    break;
	  case LIMITER_ZERO :
	    phi_ = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    phi_ = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
					  u0Min, u0Max, 4);
	    break;
	  case LIMITER_VANLEER :
	    phi_ = Limiter_VanLeer(uQuad, W[i][j][n], 
				  u0Min, u0Max, 4);
	    break;
	  case LIMITER_VANALBADA :
	    phi_ = Limiter_VanAlbada(uQuad, W[i][j][n], 
				    u0Min, u0Max, 4);
	    break;
	  default:
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  } /* endswitch */
	    
	  phi[i][j][n] = phi_;
	} /* endfor */
      } // end limiter if
	  
    } else {
      dWdx[i][j].Vacuum();
      dWdy[i][j].Vacuum();
      phi[i][j].Vacuum();
    } /* endif */

}

/********************************************************
 * Routine: Linear_Reconstruction_GreenGauss            *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh for the specified quadrilateral   *
 * solution block.  A Green-Gauss approach is used      *
 * in the evaluation of the unlimited solution          *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Flame2D_Quad_Block::Linear_Reconstruction_GreenGauss(const int Limiter) {

  // Compute and store the primitive state at the nodes.
  Update_Nodal_Values();

  // Carry out the limited solution reconstruction in
  // each cell of the computational mesh.
  for (int j  = JCl-Nghost+1 ; j <= JCu+Nghost-1 ; ++j ) {
    for ( int i = ICl-Nghost+1 ; i <= ICu+Nghost-1 ; ++i ) {
      Linear_Reconstruction_GreenGauss(i, j, Limiter);
    } /* endfor */
  } /* endfor */

}



void Flame2D_Quad_Block::Linear_Reconstruction_LeastSquares(const int i, 
							    const int j,
							    const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi_;
    double DxDx_ave, DxDy_ave, DyDy_ave, DU;
    static Vector2D dX, dXe, dXw, dXn, dXs;
    static Flame2D_State DUDx_ave, DUDy_ave;
    
    const int NUM_VAR_FLAME2D(NumVar());

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */
    
    /****************************************************************/
    //FOR VISCOUS -> CHANGED TO USE ALL 8
    if (i == ICl-Nghost || i == ICu+Nghost ||
	j == JCl-Nghost || j == JCu+Nghost) {
      n_pts = 0;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j-1;
      i_index[3] = i-1; j_index[3] = j  ;
      i_index[4] = i+1; j_index[4] = j  ;
      i_index[5] = i-1; j_index[5] = j+1;
      i_index[6] = i  ; j_index[6] = j+1;
      i_index[7] = i+1; j_index[7] = j+1;
    }    
    /****************************************************************/
    
    if (n_pts > 0) {
      DUDx_ave.Vacuum();
      DUDy_ave.Vacuum();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
     
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc;
	dX -= Grid.Cell[i][j].Xc;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
	for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	  DU = W[ i_index[n2] ][ j_index[n2] ][k] - W[i][j][k];
	  DUDx_ave[k] += DU*dX.x;
	  DUDy_ave[k] += DU*dX.y;
	}
      } /* endfor */
      
      // don't need to do this, it will cancel out
      // DUDx_ave /= double(n_pts);
      // DUDy_ave /= double(n_pts);
      // DxDx_ave /= double(n_pts);
      // DxDy_ave /= double(n_pts);
      // DyDy_ave /= double(n_pts);
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	dWdx[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
	dWdy[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      }

      if (!Freeze_Limiter) {

	dXe = Grid.xfaceE(i, j) - Grid.Cell[i][j].Xc;
	dXw = Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc;
	dXn = Grid.xfaceN(i, j) - Grid.Cell[i][j].Xc;
	dXs = Grid.xfaceS(i, j) - Grid.Cell[i][j].Xc;

	for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
	  u0Min = W[i][j][n];
	  u0Max = u0Min;
	  for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	    u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	    u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
	  } /* endfor */
	  
	  uQuad[0] = W[i][j][n] + 
	    dWdx[i][j][n]*dXe.x +
	    dWdy[i][j][n]*dXe.y ;
	  uQuad[1] = W[i][j][n] + 
	    dWdx[i][j][n]*dXw.x +
	    dWdy[i][j][n]*dXw.y ;
	  uQuad[2] = W[i][j][n] + 
	    dWdx[i][j][n]*dXn.x +
	    dWdy[i][j][n]*dXn.y ;
	  uQuad[3] = W[i][j][n] + 
	    dWdx[i][j][n]*dXs.x +
	    dWdy[i][j][n]*dXs.y ;
	  
	  switch(Limiter) {
	  case LIMITER_ONE :
	    phi_ = ONE;
	    break;
	  case LIMITER_ZERO :
	    phi_ = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    phi_ = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
					  u0Min, u0Max, 4);
	    break;
	  case LIMITER_VANLEER :
	    phi_ = Limiter_VanLeer(uQuad, W[i][j][n], 
				  u0Min, u0Max, 4);
	    break;
	  case LIMITER_VANALBADA :
	    phi_ = Limiter_VanAlbada(uQuad, W[i][j][n], 
				    u0Min, u0Max, 4);
	    break;
	  default:
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  } /* endswitch */
	    
	  phi[i][j][n] = phi_;
	} /* endfor */
      }//end limiter if
    } else {
        dWdx[i][j].Vacuum();
        dWdy[i][j].Vacuum();
        phi[i][j].Vacuum();  
    } /* endif */
    
 

}


/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh of the specified quadrilateral    *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Flame2D_Quad_Block::Linear_Reconstruction_LeastSquares(const int Limiter) {
 
    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

    for (int j  = JCl-Nghost+1 ; j <= JCu+Nghost-1 ; ++j ) {
       for (int i = ICl-Nghost+1 ; i <= ICu+Nghost-1 ; ++i ) {
	   Linear_Reconstruction_LeastSquares(i, j, Limiter);
           //Linear_Reconstruction_LeastSquares_2(SolnBlk, i, j, Limiter);
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares_2        *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within a given cell (i,j) of   *
 * the computational mesh for the specified             *
 * quadrilateral solution block.  A least squares       *
 * approach is used in the evaluation of the unlimited  *
 * solution gradients.  Several slope limiters may be   *
 * used.                                                *
 *                                                      *
 ********************************************************/
void Flame2D_Quad_Block::Linear_Reconstruction_LeastSquares_2(const int i, 
							      const int j,
							      const int Limiter) {

    int n, n2, n_pts, i_index[8], j_index[8];
    double u0Min, u0Max, uQuad[4], phi_;
    double DU, DxDx_ave, DxDy_ave, DyDy_ave;
    static Vector2D dX, dXe, dXw, dXn, dXs;
    static Flame2D_State DUDx_ave, DUDy_ave;
   
    const int NUM_VAR_FLAME2D(NumVar());

    /* Carry out the limited solution reconstruction in
       each cell of the computational mesh. */

     if (i == ICl-Nghost || i == ICu+Nghost ||
	 j == JCl-Nghost || j == JCu+Nghost) {
       n_pts = 0;
    } else if ((i == ICl-Nghost+1) && 
               (Grid.BCtypeW[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
                 Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
                 Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
                 Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
		 Grid.BCtypeW[j] == BC_1DFLAME_OUTFLOW) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (j == JCl) {
            n_pts = 0;
         } else if (j == JCu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */           
    } else if ((i == ICu+Nghost-1) && 
               (Grid.BCtypeE[j] != BC_NONE)) {
      if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
         n_pts = 0;
      } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
                 Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
                 Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
                 Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		 Grid.BCtypeE[j] == BC_1DFLAME_OUTFLOW) {
         if (j == JCl) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j  ;
            i_index[1] = i+1; j_index[1] = j  ;
            i_index[2] = i-1; j_index[2] = j+1;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (j == JCu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (j == JCl) {
            n_pts = 0;
         } else if (j == JCu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
    } else if ((j == JCl-Nghost+1) && 
               (Grid.BCtypeS[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
                 Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
                 Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
                 Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
		 Grid.BCtypeS[i] == BC_1DFLAME_OUTFLOW) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (i == ICl) {
            n_pts = 0;
         } else if (i == ICu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
    } else if ((j == JCu+Nghost-1) && 
               (Grid.BCtypeN[i] != BC_NONE)) {
      if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
         n_pts = 0;
      } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
                 Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
                 Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
                 Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		 Grid.BCtypeN[i] == BC_1DFLAME_OUTFLOW) {
         if (i == ICl) {
            n_pts = 5;
            i_index[0] = i  ; j_index[0] = j-1;
            i_index[1] = i+1; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j  ;
            i_index[3] = i  ; j_index[3] = j+1;
            i_index[4] = i+1; j_index[4] = j+1;
         } else if (i == ICu) {
            n_pts = 5;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i-1; j_index[2] = j  ;
            i_index[3] = i-1; j_index[3] = j+1;
            i_index[4] = i  ; j_index[4] = j+1;
         } else {
            n_pts = 8;
            i_index[0] = i-1; j_index[0] = j-1;
            i_index[1] = i  ; j_index[1] = j-1;
            i_index[2] = i+1; j_index[2] = j-1;
            i_index[3] = i-1; j_index[3] = j  ;
            i_index[4] = i+1; j_index[4] = j  ;
            i_index[5] = i-1; j_index[5] = j+1;
            i_index[6] = i  ; j_index[6] = j+1;
            i_index[7] = i+1; j_index[7] = j+1;
         } /* endif */
      } else {
         if (i == ICl) {
            n_pts = 0;
         } else if (i == ICu) {
            n_pts = 0;
         } else {
            n_pts = 0;
         } /* endif */
      } /* endif */
    } else if ((i == ICl && j == JCl) && 
               (Grid.BCtypeW[j] != BC_NONE &&
                Grid.BCtypeS[i] != BC_NONE)) {
      n_pts = 3;
      i_index[0] = i+1; j_index[0] = j  ;
      i_index[1] = i  ; j_index[1] = j+1;
      i_index[2] = i+1; j_index[2] = j+1;
    } else if ((i == ICl && j == JCu) && 
               (Grid.BCtypeW[j] != BC_NONE &&
                Grid.BCtypeN[i] != BC_NONE)) {
      n_pts = 3;
      i_index[0] = i  ; j_index[0] = j-1;
      i_index[1] = i+1; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j  ;
    } else if ((i == ICu && j == JCl) && 
               (Grid.BCtypeE[j] != BC_NONE &&
                Grid.BCtypeS[i] != BC_NONE)) {
      n_pts = 3;
      i_index[0] = i-1; j_index[0] = j  ;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i  ; j_index[2] = j+1;
    } else if ((i == ICu && j == JCu) && 
               (Grid.BCtypeE[j] != BC_NONE &&
                Grid.BCtypeN[i] != BC_NONE)) {
      n_pts = 3;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i-1; j_index[2] = j  ;
    } else if ((i == ICl) && 
               (Grid.BCtypeW[j] != BC_NONE)) {
      n_pts = 5;
      i_index[0] = i  ; j_index[0] = j-1;
      i_index[1] = i+1; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j  ;
      i_index[3] = i  ; j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
    } else if ((j == JCu) && 
               (Grid.BCtypeN[i] != BC_NONE)) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j-1;
      i_index[3] = i-1; j_index[3] = j  ;
      i_index[4] = i+1; j_index[4] = j  ;
    } else if ((i == ICu) && 
               (Grid.BCtypeE[j] != BC_NONE)) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i-1; j_index[2] = j  ;
      i_index[3] = i-1; j_index[3] = j+1;
      i_index[4] = i  ; j_index[4] = j+1;

    } else if ((j == JCl) && 
               (Grid.BCtypeS[i] != BC_NONE)) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j  ;
      i_index[1] = i+1; j_index[1] = j  ;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i  ; j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i  ; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j-1;
      i_index[3] = i-1; j_index[3] = j  ;
      i_index[4] = i+1; j_index[4] = j  ;
      i_index[5] = i-1; j_index[5] = j+1;
      i_index[6] = i  ; j_index[6] = j+1;
      i_index[7] = i+1; j_index[7] = j+1;
    } /* endif */
    
    if (n_pts > 0) {
        DUDx_ave.Vacuum();
        DUDy_ave.Vacuum();
        DxDx_ave = ZERO;
        DxDy_ave = ZERO;
        DyDy_ave = ZERO;
    
        for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	  dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc;
	  dX -= Grid.Cell[i][j].Xc;
	  DxDx_ave += dX.x*dX.x;
	  DxDy_ave += dX.x*dX.y;
	  DyDy_ave += dX.y*dX.y;
	  for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	    DU = W[ i_index[n2] ][ j_index[n2] ][k] - W[i][j][k];
	    DUDx_ave[k] += DU*dX.x;
	    DUDy_ave[k] += DU*dX.y;
	  }
        } /* endfor */
    					    
      // don't need to do this, it will cancel out
      // DUDx_ave /= double(n_pts);
      // DUDy_ave /= double(n_pts);
      // DxDx_ave /= double(n_pts);
      // DxDy_ave /= double(n_pts);
      // DyDy_ave /= double(n_pts);
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	dWdx[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
	dWdy[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				  (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      }

      // Calculate slope limiters. 
      if (!Freeze_Limiter) {

	dXe = Grid.xfaceE(i, j) - Grid.Cell[i][j].Xc;
	dXw = Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc;
	dXn = Grid.xfaceN(i, j) - Grid.Cell[i][j].Xc;
	dXs = Grid.xfaceS(i, j) - Grid.Cell[i][j].Xc;

	for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
	  u0Min = W[i][j][n];
	  u0Max = u0Min;
	  for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	    u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	    u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
	  } /* endfor */
	    
	  uQuad[0] = W[i][j][n] + 
	    dWdx[i][j][n]*dXe.x +
	    dWdy[i][j][n]*dXe.y ;
	  uQuad[1] = W[i][j][n] + 
	    dWdx[i][j][n]*dXw.x +
	    dWdy[i][j][n]*dXw.y ;
	  uQuad[2] = W[i][j][n] + 
	    dWdx[i][j][n]*dXn.x +
	    dWdy[i][j][n]*dXn.y ;
	  uQuad[3] = W[i][j][n] + 
	    dWdx[i][j][n]*dXs.x +
	    dWdy[i][j][n]*dXs.y ;
	    
	  switch(Limiter) {
	  case LIMITER_ONE :
	    phi_ = ONE;
	    break;
	  case LIMITER_ZERO :
	    phi_ = ZERO;
	    break;
	  case LIMITER_BARTH_JESPERSEN :
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  case LIMITER_VENKATAKRISHNAN :
	    phi_ = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
					  u0Min, u0Max, 4);	  
	    break;
	  case LIMITER_VANLEER :
	    phi_ = Limiter_VanLeer(uQuad, W[i][j][n], 
				  u0Min, u0Max, 4);
	    break;
	  case LIMITER_VANALBADA :
	    phi_ = Limiter_VanAlbada(uQuad, W[i][j][n], 
				    u0Min, u0Max, 4);
	    break;
	  default:
	    phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					 u0Min, u0Max, 4);
	    break;
	  } /* endswitch */
	    
	  phi[i][j][n] = phi_;
	} /* endfor */
      } //end limiter if

    } else {
      dWdx[i][j].Vacuum();
      dWdy[i][j].Vacuum();
      phi[i][j].Vacuum();
    } /* endif */
   
}

/////////////////////////////////////////////////////////////////////
/// Viscous Reconstruction Methods
/////////////////////////////////////////////////////////////////////

//Diamond Path
/********************************************************
 * Routine: Linear_Reconstruction_LeastSquares          *
 *                                                      *
 * Performs the reconstruction of a limited piecewise   *
 * linear solution state within each cell of the        *
 * computational mesh of the specified quadrilateral    *
 * solution block.  A least squares approach is         *
 * used in the evaluation of the unlimited solution     *
 * gradients.  Several slope limiters may be used.      *
 *                                                      *
 ********************************************************/
void Flame2D_Quad_Block::Linear_Reconstruction_LeastSquares_Diamond(const int Limiter) {
 

  // Compute and store the primitive state at the nodes.
  Update_Nodal_Values();

  // Carry out the limited solution reconstruction in
  // each cell of the computational mesh.
  
  for (int j  = JCl-Nghost+1 ; j <= JCu+Nghost-1 ; ++j ) {
    for (int i = ICl-Nghost+1 ; i <= ICu+Nghost-1 ; ++i ) {
      Linear_Reconstruction_LeastSquares_Diamond(i, j, Limiter);
      
    } /* endfor */
  } /* endfor */
  
}


//Diamond path 
void Flame2D_Quad_Block::Linear_Reconstruction_LeastSquares_Diamond(const int i, 
								    const int j,
								    const int Limiter) {


  int n, n2, n_pts, i_index[8], j_index[8];
  int n_neigbour;
  double u0Min, u0Max, uQuad[4], phi_;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double area[4];

  static Vector2D dXe, dXw, dXn, dXs, dX[4];  
  static Flame2D_State DU[4], DUDx_ave, DUDy_ave;
  double QuadraturePoint_N, 
    QuadraturePoint_E, 
    QuadraturePoint_S, 
    QuadraturePoint_W;  
  const Flame2D_State *TopVertex, *BottomVertex;

  const int NUM_VAR_FLAME2D = NumVar();

  /****************************************************************/
  //  * A least squares       *
  //  * approach is used in the evaluation of the unlimited  *
  //  * solution gradients on cell faces that is peformed on a diamond path
  // For cell (i,j), in order to reconstruct the gradients on cell faces, 
  // information of four points are needed
  // centeroid (i+1, j) and vertices (i+1, j) and (i+1, j+1)
  // By convention: the quadrature point (mid-edge) (i, j)  

  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
 
 
  } else {
 
     n_pts = 8;
     i_index[0] = i-1; j_index[0] = j-1;
     i_index[1] = i  ; j_index[1] = j-1;
     i_index[2] = i+1; j_index[2] = j-1;
     i_index[3] = i-1; j_index[3] = j  ;
     i_index[4] = i+1; j_index[4] = j  ;
     i_index[5] = i-1; j_index[5] = j+1;
     i_index[6] = i  ; j_index[6] = j+1;
     i_index[7] = i+1; j_index[7] = j+1;

  }    

  n_neigbour = 4;
  
  if (n_pts > 0) {
 
    /*************** NORTH ****************************/
    /*Formulate the gradients of primitive parameters on the north face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
 

    TopVertex = &Wnd[i][j+1];
    BottomVertex = &Wnd[i+1][j+1];

    dX[0] = Grid.Cell[i][j].Xc - Grid.xfaceN(i,j);
    dX[1] = Grid.Cell[i][j+1].Xc - Grid.xfaceN(i,j);
    dX[2] = Grid.Node[i][j+1].X -  Grid.xfaceN(i,j);
    dX[3] = Grid.Node[i+1][j+1].X -  Grid.xfaceN(i,j);
    //The calculation of this area is used to weight the gradient at the cell center      
    area[0] = HALF*(Grid.Node[i+1][j+1].X - Grid.Node[i][j+1].X )^
      (Grid.xfaceN(i,j)- Grid.Cell[i][j].Xc);
  
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      QuadraturePoint_N = HALF*( (*TopVertex)[k] + (*BottomVertex)[k]);
      //Left state cell (i,j)
      DU[0][k] = W[i][j][k] - QuadraturePoint_N;
      //rigth state cell (i, j+1)
      DU[1][k] = W[i][j+1][k] - QuadraturePoint_N;
      //top vertex  (i, j+1)
      DU[2][k] = (*TopVertex)[k] - QuadraturePoint_N;
      //bottom vertex (i+1,j+1) 
      DU[3][k] = (*BottomVertex)[k] - QuadraturePoint_N;
    }

 
    for ( n2 = 0 ; n2 <= n_neigbour -1 ; ++n2 ) {
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	DUDx_ave[k] += DU[n2][k]*dX[n2].x;
	DUDy_ave[k] += DU[n2][k]*dX[n2].y;
      }
    } /* endfor */
 
    // don't need to do this, it will cancel out
    // DUDx_ave = DUDx_ave/double(n_neigbour);
    // DUDy_ave = DUDy_ave/double(n_neigbour);
    // DxDx_ave = DxDx_ave/double(n_neigbour);
    // DxDy_ave = DxDy_ave/double(n_neigbour);
    // DyDy_ave = DyDy_ave/double(n_neigbour);
 
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx_faceN[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      dWdy_faceN[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    }
     
    /*************** EAST *****************************/
    /*Formulate the gradients of primitive parameters on the east face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
   
    //needs to assign topvertex and bottomvertex information
    TopVertex = &Wnd[i+1][j+1];
    BottomVertex =  &Wnd[i+1][j];

    dX[0] = Grid.Cell[i][j].Xc - Grid.xfaceE(i, j);
    dX[1] = Grid.Cell[i+1][j].Xc - Grid.xfaceE(i,j);
    dX[2] = Grid.Node[i+1][j+1].X - Grid.xfaceE(i,j);
    dX[3] = Grid.Node[i+1][j].X - Grid.xfaceE(i,j);
    //The calculation of this area is used to weight the gradient at the cell center      
    area[1] = HALF*(Grid.Node[i+1][j+1].X - Grid.Node[i+1][j].X )^
      (Grid.Cell[i][j].Xc - Grid.xfaceE(i,j));

    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      QuadraturePoint_E = HALF*( (*TopVertex)[k] + (*BottomVertex)[k]);
      //Left state cell (i,j)
      DU[0][k] = W[i][ j][k] - QuadraturePoint_E;
      //rigth state cell (i+1, j)
      DU[1][k] = W[i+1][ j][k] - QuadraturePoint_E;
      //top vertex  (i+1, j+1)
      DU[2][k] = (*TopVertex)[k] - QuadraturePoint_E;
      //bottom vertex (i+1,j) 
      DU[3][k] = (*BottomVertex)[k] - QuadraturePoint_E;
    }
 
    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	DUDx_ave[k] += DU[n2][k]*dX[n2].x;
	DUDy_ave[k] += DU[n2][k]*dX[n2].y;
      }      
    } /* endfor */
 
    // don't need to do this, it will cancel out
    // DUDx_ave = DUDx_ave/double(n_neigbour);
    // DUDy_ave = DUDy_ave/double(n_neigbour);
    // DxDx_ave = DxDx_ave/double(n_neigbour);
    // DxDy_ave = DxDy_ave/double(n_neigbour);
    // DyDy_ave = DyDy_ave/double(n_neigbour);
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx_faceE[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      dWdy_faceE[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    }
     
    /*************** WEST *****************************/
    /*Formulate the gradients of primitive parameters on the west face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    TopVertex = &Wnd[i][j];
    BottomVertex =  &Wnd[i][j+1];

    dX[0] = Grid.Cell[i][j].Xc - Grid.xfaceW(i,j);
    dX[1] = Grid.Cell[i-1][j].Xc - Grid.xfaceW(i,j);
    dX[2] = Grid.Node[i][j].X - Grid.xfaceW(i,j);
    dX[3] = Grid.Node[i][j+1].X - Grid.xfaceW(i,j);
    //The calculation of this area is used to weight the gradient at the cell center      
    area[2] = HALF*(Grid.Node[i][j+1].X - Grid.Node[i][j].X )^
      ( Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc );

    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      QuadraturePoint_W = HALF*( (*TopVertex)[k] + (*BottomVertex)[k]);
      //Left state cell (i,j)
      DU[0][k] = W[i][ j][k] - QuadraturePoint_W;
      //rigth state cell (i-1, j)
      DU[1][k] = W[i-1][j][k] - QuadraturePoint_W;
      //top vertex  (i, j)
      DU[2][k] = (*TopVertex)[k] - QuadraturePoint_W;
      //bottom vertex (i,j+1) 
      DU[3][k] = (*BottomVertex)[k] - QuadraturePoint_W;
    }

    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {      
	DUDx_ave[k] += DU[n2][k]*dX[n2].x;
	DUDy_ave[k] += DU[n2][k]*dX[n2].y;
      }
    } /* endfor */
 
     // don't need to do this, it will cancel out
    // DUDx_ave = DUDx_ave/double(n_neigbour);
    // DUDy_ave = DUDy_ave/double(n_neigbour);
    // DxDx_ave = DxDx_ave/double(n_neigbour);
    // DxDy_ave = DxDy_ave/double(n_neigbour);
    // DyDy_ave = DyDy_ave/double(n_neigbour);
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx_faceW[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      dWdy_faceW[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    }

    /*************** SOUTH ****************************/
    /*Formulate the gradients of primitive parameters on the south face of cell (i, j)*/
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    TopVertex = &Wnd[i+1][j];
    BottomVertex =  &Wnd[i][j];

    dX[0] = Grid.Cell[i][j].Xc - Grid.xfaceS(i, j);
    dX[1] = Grid.Cell[i][j-1].Xc - Grid.xfaceS(i, j);
    dX[2] = Grid.Node[i+1][j].X - Grid.xfaceS(i, j);
    dX[3] = Grid.Node[i][j].X - Grid.xfaceS(i, j);
    //The calculation of this area is used to weight the gradient at the cell center      
    area[3] = HALF*(Grid.Node[i+1][j].X - Grid.Node[i][j].X )^
      (Grid.Cell[i][j].Xc - Grid.xfaceS(i,j) );

    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      QuadraturePoint_S = HALF*((*TopVertex)[k] + (*BottomVertex)[k]);
      //Left state cell (i,j)
      DU[0][k] = W[i][j][k] - QuadraturePoint_S;
      //rigth state cell (i, j-1)
      DU[1][k] = W[i][j-1][k] - QuadraturePoint_S;
      //top vertex  (i+1, j)
      DU[2][k] = (*TopVertex)[k] - QuadraturePoint_S;
      //bottom vertex (i,j) 
      DU[3][k] = (*BottomVertex)[k] - QuadraturePoint_S;
    }

    for ( n2 = 0 ; n2 <=  n_neigbour-1 ; ++n2 ) {
      DxDx_ave += dX[n2].x*dX[n2].x;
      DxDy_ave += dX[n2].x*dX[n2].y;
      DyDy_ave += dX[n2].y*dX[n2].y;
      for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
	DUDx_ave[k] += DU[n2][k]*dX[n2].x;
	DUDy_ave[k] += DU[n2][k]*dX[n2].y;
      }      
    } /* endfor */
 
    // don't need to do this, it will cancel out
    // DUDx_ave = DUDx_ave/double(n_neigbour);
    // DUDy_ave = DUDy_ave/double(n_neigbour);
    // DxDx_ave = DxDx_ave/double(n_neigbour);
    // DxDy_ave = DxDy_ave/double(n_neigbour);
    // DyDy_ave = DyDy_ave/double(n_neigbour);
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx_faceS[i][j][k] = ( (DUDx_ave[k]*DyDy_ave-DUDy_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      dWdy_faceS[i][j][k] = ( (DUDy_ave[k]*DxDx_ave-DUDx_ave[k]*DxDy_ave)/
				      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    }

    /**************************************************/
    //Area weighted gradients at cell centers
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx[i][j][k] = ( dWdx_faceN[i][j][k]*area[0] + 
				dWdx_faceE[i][j][k]*area[1] +
				dWdx_faceW[i][j][k]*area[2] + 
				dWdx_faceS[i][j][k]*area[3] ); 
      dWdx[i][j][k] /= Grid.Cell[i][j].A; 
   
      dWdy[i][j][k] = ( dWdy_faceN[i][j][k]*area[0] + 
				dWdy_faceE[i][j][k]*area[1] +
				dWdy_faceW[i][j][k]*area[2] + 
				dWdy_faceS[i][j][k]*area[3] ); 
      dWdy[i][j][k] /= Grid.Cell[i][j].A;
    }
    /**************************************************/
   
    if (!Freeze_Limiter) {

      dXe = Grid.xfaceE(i, j) - Grid.Cell[i][j].Xc;
      dXw = Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc;
      dXn = Grid.xfaceN(i, j) - Grid.Cell[i][j].Xc;
      dXs = Grid.xfaceS(i, j) - Grid.Cell[i][j].Xc;
   
      for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for(n2 = 0 ; n2 <= n_pts-1 ; ++n2){
           u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
           u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]); 
	  
	}	
	uQuad[0] = W[i][j][n] + 
	  dWdx[i][j][n]*dXe.x +
	  dWdy[i][j][n]*dXe.y ;
	uQuad[1] = W[i][j][n] + 
	  dWdx[i][j][n]*dXw.x +
	  dWdy[i][j][n]*dXw.y ;
	uQuad[2] = W[i][j][n] + 
	  dWdx[i][j][n]*dXn.x +
	  dWdy[i][j][n]*dXn.y ;
	uQuad[3] = W[i][j][n] + 
	  dWdx[i][j][n]*dXs.x +
	  dWdy[i][j][n]*dXs.y ;
	  
	switch(Limiter) {
	case LIMITER_ONE :
	  phi_ = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_ = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
				       u0Min, u0Max, 4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_ = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
					u0Min, u0Max, 4);
	  break;
	case LIMITER_VANLEER :
	  phi_ = Limiter_VanLeer(uQuad, W[i][j][n], 
				u0Min, u0Max, 4);
	  break;
	case LIMITER_VANALBADA :
	  phi_ = Limiter_VanAlbada(uQuad, W[i][j][n], 
				  u0Min, u0Max, 4);
	  break;
	default:
	  phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
				       u0Min, u0Max, 4);
	  break;
	} /* endswitch */
	
	phi[i][j][n] = phi_;
      } /* endfor */
    }//end limiter if 
  } else {
    dWdx_faceN[i][j].Vacuum();
    dWdx_faceS[i][j].Vacuum();
    dWdx_faceE[i][j].Vacuum();
    dWdx_faceW[i][j].Vacuum();
    dWdy_faceN[i][j].Vacuum();
    dWdy_faceS[i][j].Vacuum();
    dWdy_faceE[i][j].Vacuum();
    dWdy_faceW[i][j].Vacuum();
    dWdx[i][j].Vacuum();
    dWdy[i][j].Vacuum();
    phi[i][j].Vacuum(); 
  } /* endif */
  
  
    
}

void Flame2D_Quad_Block::Linear_Reconstruction_GreenGauss_Diamond(const int Limiter) { 

  // Compute and store the primitive state at the nodes.
  Update_Nodal_Values();

  // Carry out the limited solution reconstruction in
  // each cell of the computational mesh.
  for (int j  = JCl-Nghost+1 ; j <= JCu+Nghost-1 ; ++j ) {
    for (int i = ICl-Nghost+1 ; i <= ICu+Nghost-1 ; ++i ) {
      Linear_Reconstruction_GreenGauss_Diamond(i, j, Limiter);	
    } 
  } 

}

void Flame2D_Quad_Block::Linear_Reconstruction_GreenGauss_Diamond(const int i, 
								  const int j,
								  const int Limiter) {
  
  int n_pts, i_index[8], j_index[8];
  double area[4], AREA;
  double W_average[4];
  static Vector2D norm[4]; 
  const int NUM_VAR_FLAME2D = NumVar();

  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  }   
  
  if (n_pts > 0) {
    
    const Flame2D_State &W_NE = Wnd[i+1][j+1];
    const Flame2D_State &W_NW = Wnd[i][j+1];
    const Flame2D_State &W_SW = Wnd[i][j];
    const Flame2D_State &W_SE = Wnd[i+1][j];
    
    /*************** NORTH ****************************/
    
    /*Formulate the gradients of primitive parameters on the north face of cell (i, j)*/
    
    //  normal vector of the SE side of a diamond 
    norm[0].x = Grid.nodeNE(i,j).X.y - Grid.Cell[i][j].Xc.y;
    norm[0].y = Grid.Cell[i][j].Xc.x - Grid.nodeNE(i,j).X.x;
    //  normal vector of the NE side of a diamond 
    norm[1].x = Grid.Cell[i][j+1].Xc.y - Grid.nodeNE(i,j).X.y;
    norm[1].y = Grid.nodeNE(i,j).X.x - Grid.Cell[i][j+1].Xc.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x = Grid.nodeNW(i,j).X.y - Grid.Cell[i][j+1].Xc.y;
    norm[2].y = Grid.Cell[i][j+1].Xc.x - Grid.nodeNW(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = Grid.Cell[i][j].Xc.y - Grid.nodeNW(i,j).X.y ;
    norm[3].y = Grid.nodeNW(i,j).X.x - Grid.Cell[i][j].Xc.x ;
    
    AREA =  HALF*(fabs((Grid.nodeNE(i,j).X-Grid.Cell[i][j].Xc)^
		       (Grid.nodeNW(i,j).X-Grid.Cell[i][j].Xc)) +
		  fabs((Grid.nodeNW(i,j).X-Grid.Cell[i][j+1].Xc)^
		       (Grid.nodeNE(i,j).X-Grid.Cell[i][j+1].Xc)));
    
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      // counterclockwise, starting from the cell center (i,j), nodeNE (i, j), 
      // top cell center (i, j+1), NodeNW (i, j);
      W_average[0] = HALF*(W[i][j][k]   + W_NE[k]);
      W_average[1] = HALF*(W[i][j+1][k] + W_NE[k]);
      W_average[2] = HALF*(W[i][j+1][k] + W_NW[k]);
      W_average[3] = HALF*(W[i][j][k]   + W_NW[k]);
      
      dWdx_faceN[i][j][k] = ( W_average[0]* norm[0].x +
				      W_average[1]* norm[1].x + 
				      W_average[2]* norm[2].x +
				      W_average[3]* norm[3].x)/AREA;
      dWdy_faceN[i][j][k] = ( W_average[0]* norm[0].y +
				      W_average[1]* norm[1].y + 
				      W_average[2]* norm[2].y +
				      W_average[3]* norm[3].y)/AREA;  
    }
    /*************** EAST ****************************/
    
    /*Formulate the gradients of primitive parameters on the east face of cell (i, j)*/
    
    //  normal vector of the NW side of a diamond  = - SE of previous
    norm[2] = - norm[0];
    //  normal vector of the SE side of a diamond 
    norm[0].x = Grid.Cell[i+1][j].Xc.y - Grid.nodeSE(i,j).X.y;
    norm[0].y = Grid.nodeSE(i,j).X.x - Grid.Cell[i+1][j].Xc.x;
    //  normal vector of the NE side of a diamond 
    norm[1].x = Grid.nodeNE(i,j).X.y -  Grid.Cell[i+1][j].Xc.y ;
    norm[1].y = Grid.Cell[i+1][j].Xc.x - Grid.nodeNE(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = Grid.nodeSE(i,j).X.y - Grid.Cell[i][j].Xc.y;
    norm[3].y = Grid.Cell[i][j].Xc.x - Grid.nodeSE(i,j).X.x;
    
    AREA =  HALF*(fabs((Grid.nodeNE(i,j).X-Grid.Cell[i+1][j].Xc)^
		       (Grid.nodeSE(i,j).X-Grid.Cell[i+1][j].Xc)) +
		  fabs((Grid.nodeSE(i,j).X-Grid.Cell[i][j].Xc)^
		       (Grid.nodeNE(i,j).X-Grid.Cell[i][j].Xc)));
    
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      // counterclockwise, starting from  nodeSE(i,j), cell (i+1, j), 
      // nodeNE(i,j), cell center (i, j+1)     
      W_average[2] = HALF*(W[i][j][k]    + W_NE[k]);
      W_average[0] = HALF*(W[i+1][j][k]  + W_SE[k]);      
      W_average[1] = HALF*(W[i+1][j][k]  + W_NE[k]);     
      W_average[3] = HALF*(W[i][j][k]    + W_SE[k]);
      
      dWdx_faceE[i][j][k] = ( W_average[0]* norm[0].x +
				      W_average[1]* norm[1].x + 
				      W_average[2]* norm[2].x + 
				      W_average[3]* norm[3].x)/AREA;
      dWdy_faceE[i][j][k] = ( W_average[0]* norm[0].y +
				      W_average[1]* norm[1].y + 
				      W_average[2]* norm[2].y + 
				      W_average[3]* norm[3].y)/AREA;  
    }
    /*************** SOUTH ****************************/
    
    /*Formulate the gradients of primitive parameters on the south face of cell (i, j)*/
    
    //  normal vector of the NE side of a diamond = -SW of previous
    norm[1] = -norm[3];
    //  normal vector of the SE side of a diamond 
    norm[0].x = Grid.nodeSE(i,j).X.y - Grid.Cell[i][j-1].Xc.y;
    norm[0].y = Grid.Cell[i][j-1].Xc.x - Grid.nodeSE(i,j).X.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x = Grid.nodeSW(i,j).X.y - Grid.Cell[i][j].Xc.y ;
    norm[2].y = Grid.Cell[i][j].Xc.x - Grid.nodeSW(i,j).X.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x = Grid.Cell[i][j-1].Xc.y - Grid.nodeSW(i,j).X.y;
    norm[3].y = Grid.nodeSW(i,j).X.x - Grid.Cell[i][j-1].Xc.x;
    
    AREA =  HALF*(fabs((Grid.nodeSE(i,j).X-Grid.Cell[i][j-1].Xc)^
		       (Grid.nodeSW(i,j).X-Grid.Cell[i][j-1].Xc)) +
		  fabs((Grid.nodeSE(i,j).X-Grid.Cell[i][j].Xc)^
		       (Grid.nodeSW(i,j).X-Grid.Cell[i][j].Xc)));

    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      // counterclockwise, starting from  cell (i, j-1), nodeSE(i,j) 
      // cell(i,j), nodeSW(i,j)     
      W_average[1] = HALF*(W[i][j][k]    + W_SE[k]);
      W_average[0] = HALF*(W[i][j-1][k]  + W_SE[k]);
      W_average[2] = HALF*(W[i][j][k]    + W_SW[k]);
      W_average[3] = HALF*(W[i][j-1][k]  + W_SW[k]);
      
      dWdx_faceS[i][j][k] = ( W_average[0]* norm[0].x +
				      W_average[1]* norm[1].x + 
				      W_average[2]* norm[2].x + 
				      W_average[3]* norm[3].x)/AREA;
      dWdy_faceS[i][j][k] = ( W_average[0]* norm[0].y +
				      W_average[1]* norm[1].y + 
				      W_average[2]* norm[2].y + 
				      W_average[3]* norm[3].y)/AREA;  
    }
    /*************** WEST ****************************/    
    /*Formulate the gradients of primitive parameters on the west face of cell (i, j ) */
    
    //  normal vector of the SE side of a diamond = - NW of previous
    norm[0] = - norm[2];
    //  normal vector of the NE side of a diamond 
    norm[1].x =  Grid.nodeNW(i,j).X.y - Grid.Cell[i][j].Xc.y;
    norm[1].y =  Grid.Cell[i][j].Xc.x - Grid.nodeNW(i,j).X.x;
    //  normal vector of the NW side of a diamond 
    norm[2].x =  Grid.Cell[i-1][j].Xc.y - Grid.nodeNW(i,j).X.y ;
    norm[2].y =  Grid.nodeNW(i,j).X.x - Grid.Cell[i-1][j].Xc.x;
    //  normal vector of the SW side of a diamond 
    norm[3].x =  Grid.nodeSW(i,j).X.y - Grid.Cell[i-1][j].Xc.y;
    norm[3].y =  Grid.Cell[i-1][j].Xc.x - Grid.nodeSW(i,j).X.x;
    
    AREA =  HALF*(fabs((Grid.nodeNW(i,j).X-Grid.Cell[i][j].Xc)^
		       (Grid.nodeSW(i,j).X-Grid.Cell[i][j].Xc)) +
		  fabs((Grid.nodeNW(i,j).X-Grid.Cell[i-1][j].Xc)^
		       (Grid.nodeSW(i,j).X-Grid.Cell[i-1][j].Xc)));
        
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      // counterclockwise, starting from  NodeSW(i,j) 
      // cell(i,j), nodeNW(i,j), cell (i-1, j)     
      W_average[0] = HALF*(W[i][j][k]    + W_SW[k]);
      W_average[1] = HALF*(W[i][j][k]    + W_NW[k]);
      W_average[2] = HALF*(W[i-1][j][k]  + W_NW[k]);
      W_average[3] = HALF*(W[i-1][j][k]  + W_SW[k]);
      
      dWdx_faceW[i][j][k] = ( W_average[0]* norm[0].x +
				      W_average[1]* norm[1].x + 
				      W_average[2]* norm[2].x + 
				      W_average[3]* norm[3].x)/AREA;
      dWdy_faceW[i][j][k] = ( W_average[0]* norm[0].y +
				      W_average[1]* norm[1].y + 
				      W_average[2]* norm[2].y + 
				      W_average[3]* norm[3].y)/AREA;  
    }
    /*************** CENTER ****************************/
    
    // area weighted gradients at cell centers, 4 inside triangles
    area[0] = HALF*(Grid.Node[i+1][j+1].X - Grid.Node[i][j+1].X )^
      (Grid.xfaceN(i,j)- Grid.Cell[i][j].Xc);
    area[1] = HALF*(Grid.Node[i+1][j+1].X - Grid.Node[i+1][j].X )^
      (Grid.Cell[i][j].Xc - Grid.xfaceE(i,j));
    area[2] = HALF*(Grid.Node[i][j+1].X - Grid.Node[i][j].X )^
      ( Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc );
    area[3] = HALF*(Grid.Node[i+1][j].X - Grid.Node[i][j].X )^
      (Grid.Cell[i][j].Xc - Grid.xfaceS(i,j) );
        
    //Reconstructed cell center gradients
    for (int k=1; k<=NUM_VAR_FLAME2D; k++) {
      dWdx[i][j][k] = ( dWdx_faceN[i][j][k]*area[0] + 
				dWdx_faceE[i][j][k]*area[1] +
				dWdx_faceW[i][j][k]*area[2] + 
				dWdx_faceS[i][j][k]*area[3] );
      dWdx[i][j][k] /= Grid.Cell[i][j].A; 
      
      dWdy[i][j][k] = ( dWdy_faceN[i][j][k]*area[0] + 
				dWdy_faceE[i][j][k]*area[1] +
				dWdy_faceW[i][j][k]*area[2] + 
				dWdy_faceS[i][j][k]*area[3] );
      dWdy[i][j][k] /= Grid.Cell[i][j].A; 
    }
    /****************************************************/
        
    double u0Min, u0Max, uQuad[4], phi_;
    double DxDx_ave, DxDy_ave, DyDy_ave;
    static Vector2D dXe, dXw, dXn, dXs;
    
    // Calculate slope limiters.  
    if (!Freeze_Limiter) {

      dXe = Grid.xfaceE(i, j) - Grid.Cell[i][j].Xc;
      dXw = Grid.xfaceW(i, j) - Grid.Cell[i][j].Xc;
      dXn = Grid.xfaceN(i, j) - Grid.Cell[i][j].Xc;
      dXs = Grid.xfaceS(i, j) - Grid.Cell[i][j].Xc;

      for ( int n=1 ; n <= NumVar(); n++ ) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for(int n2 = 0 ; n2 <= n_pts-1 ; ++n2){
	  u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	  u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]); 
	}
	
	uQuad[0] = W[i][j][n] + 
	  dWdx[i][j][n]*dXe.x +
	  dWdy[i][j][n]*dXe.y ;
	uQuad[1] = W[i][j][n] + 
	  dWdx[i][j][n]*dXw.x +
	  dWdy[i][j][n]*dXw.y ;
	uQuad[2] = W[i][j][n] + 
	  dWdx[i][j][n]*dXn.x +
	  dWdy[i][j][n]*dXn.y ;
	uQuad[3] = W[i][j][n] + 
	  dWdx[i][j][n]*dXs.x +
	  dWdy[i][j][n]*dXs.y ;
	  
	switch(Limiter) {
	case LIMITER_ONE :
	  phi_ = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_ = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
				       u0Min, u0Max, 4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_ = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
					u0Min, u0Max, 4);
	  break;
	case LIMITER_VANLEER :
	  phi_ = Limiter_VanLeer(uQuad, W[i][j][n], 
				u0Min, u0Max, 4);
	  break;
	case LIMITER_VANALBADA :
	  phi_ = Limiter_VanAlbada(uQuad, W[i][j][n], 
				  u0Min, u0Max, 4);
	  break;
	default:
	  phi_ = Limiter_BarthJespersen(uQuad, W[i][j][n], 
				       u0Min, u0Max, 4);
	  break;
	} /* endswitch */
	
	phi[i][j][n] = phi_;
      } /* endfor */
    }//end limiter if 
  } else {
    dWdx_faceN[i][j].Vacuum();
    dWdx_faceS[i][j].Vacuum();
    dWdx_faceE[i][j].Vacuum();
    dWdx_faceW[i][j].Vacuum();
    dWdy_faceN[i][j].Vacuum();
    dWdy_faceS[i][j].Vacuum();
    dWdy_faceE[i][j].Vacuum();
    dWdy_faceW[i][j].Vacuum();
    dWdx[i][j].Vacuum();
    dWdy[i][j].Vacuum();
    phi[i][j].Vacuum(); 
  } 
      
}

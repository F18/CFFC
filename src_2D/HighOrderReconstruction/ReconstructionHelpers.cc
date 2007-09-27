/*!\file ReconstructionHelpers.cc
  \brief Source file providing implementation of subroutines prototyped in ReconstructionHelpers.h */

#include "ReconstructionHelpers.h"

/******************************************************************************************
Function for generating the geometric coefficients for a cartesian cell
******************************************************************************************/
double GeomCoeffCartesian(int p1, int p2, double deltaX, double deltaY, double deltaXC, double deltaYC){

  /* p1 -> the first power coefficient
     p2 -> the second power coefficient
     deltaX -> the grid size in the X direction
     deltaY -> the grid size in the Y direction
     deltaXC -> the X distance between the center of the reconstructed cell and that of the cell used in the reconstruction
     deltaYC -> the Y distance between the center of the reconstructed cell and that of the cell used in the reconstruction

     Obs. To compute the coefficient of the reconstructed cell, deltaXC and deltaYC must be ZERO. 
  */

  double val1, val2;
  double coef_x1, coef_x2, coef_y1, coef_y2;

  val1 = val2 = 0.0;
  coef_x1 = deltaX/2  + deltaXC;
  coef_x2 = -deltaX/2 + deltaXC;
  coef_y1 = deltaY/2  + deltaYC;
  coef_y2 = -deltaY/2 + deltaYC;

  for (int m=1; m<=p1+1; m++){
    val1 += pow(coef_x1,p1+1-m)*pow(coef_x2,m-1);
  }
  for (int l=1; l<=p2+1; l++){
    val2 += pow(coef_y1,p2+1-l)*pow(coef_y2,l-1);
  }

  return val1*val2/((p1+1)*(p2+1));
}

// MakeReconstructionStencil(int,int,vector<int>) function
// Set the stencil for the 1D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, vector<int> & i_index){

  // Obs. The first position (i_index[0]) corresponds to (iCell)

  switch(rings){

  case 3: // three rings of cells around (iCell)
    /* Third ring */
    i_index[5]=iCell-3;
    i_index[6]=iCell+3;

  case 2: // two rings of cells around (iCell)
    /* Second ring */
    i_index[3]=iCell-2;
    i_index[4]=iCell+2;

  case 1: // one ring of cells around (iCell,jCell)
    /* First ring */
    i_index[1]=iCell-1;
    i_index[2]=iCell+1;

  case 0: 
    i_index[0]=iCell; /* cell (iCell) */

  default: // general expression
    i_index[0]=iCell; /* cell (iCell) */
    for (int i=iCell-rings, Pos=1; i<=iCell+rings; ++i){
      if(i!=iCell){
	i_index[Pos] = i;
	++Pos;
      }
    }
  }
}

// MakeReconstructionStencil(int,int,int,vector<int>,vector<int>) function
// Set the stencil for the 2D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       vector<int> & i_index, vector<int> & j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    i_index[9] =iCell-2;  j_index[9]=jCell-2;
    i_index[10]=iCell-1; j_index[10]=jCell-2;
    i_index[11]=iCell  ; j_index[11]=jCell-2;
    i_index[12]=iCell+1; j_index[12]=jCell-2;
    i_index[13]=iCell+2; j_index[13]=jCell-2;
    i_index[14]=iCell-2; j_index[14]=jCell-1;
    i_index[15]=iCell+2; j_index[15]=jCell-1;
    i_index[16]=iCell-2; j_index[16]=jCell;
    i_index[17]=iCell+2; j_index[17]=jCell;
    i_index[18]=iCell-2; j_index[18]=jCell+1;
    i_index[19]=iCell+2; j_index[19]=jCell+1;
    i_index[20]=iCell-2; j_index[20]=jCell+2;
    i_index[21]=iCell-1; j_index[21]=jCell+2;
    i_index[22]=iCell  ; j_index[22]=jCell+2;
    i_index[23]=iCell+1; j_index[23]=jCell+2;
    i_index[24]=iCell+2; j_index[24]=jCell+2;

  case 1: // one ring of cells around (iCell,jCell)

    i_index[0]=iCell;   j_index[0]=jCell; /* cell (iCell,jCell) */
    /* First ring */
    i_index[1]=iCell-1; j_index[1]=jCell-1;
    i_index[2]=iCell;   j_index[2]=jCell-1;
    i_index[3]=iCell+1; j_index[3]=jCell-1;
    i_index[4]=iCell-1; j_index[4]=jCell;
    i_index[5]=iCell+1; j_index[5]=jCell;
    i_index[6]=iCell-1; j_index[6]=jCell+1;
    i_index[7]=iCell;   j_index[7]=jCell+1;
    i_index[8]=iCell+1; j_index[8]=jCell+1;
    break;

  default: // general expression
    i_index[0] = iCell;
    j_index[0] = jCell;
    for (int i=iCell-rings, Poz=1; i<=iCell+rings; ++i)
      for (int j=jCell-rings; j<=jCell+rings; ++j){
	if(!((i==iCell)&&(j==jCell)) ){
	  i_index[Poz] = i;
	  j_index[Poz] = j;
	  ++Poz;
	}
      }
  }//endswitch
}

//MakeReconstructionStencil(int,int,int,int *,int *) function
// Set the stencil for the 2D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       int *i_index, int *j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    i_index[9] =iCell-2;  j_index[9]=jCell-2;
    i_index[10]=iCell-1; j_index[10]=jCell-2;
    i_index[11]=iCell  ; j_index[11]=jCell-2;
    i_index[12]=iCell+1; j_index[12]=jCell-2;
    i_index[13]=iCell+2; j_index[13]=jCell-2;
    i_index[14]=iCell-2; j_index[14]=jCell-1;
    i_index[15]=iCell+2; j_index[15]=jCell-1;
    i_index[16]=iCell-2; j_index[16]=jCell;
    i_index[17]=iCell+2; j_index[17]=jCell;
    i_index[18]=iCell-2; j_index[18]=jCell+1;
    i_index[19]=iCell+2; j_index[19]=jCell+1;
    i_index[20]=iCell-2; j_index[20]=jCell+2;
    i_index[21]=iCell-1; j_index[21]=jCell+2;
    i_index[22]=iCell  ; j_index[22]=jCell+2;
    i_index[23]=iCell+1; j_index[23]=jCell+2;
    i_index[24]=iCell+2; j_index[24]=jCell+2;

  case 1: // one ring of cells around (iCell,jCell)

    i_index[0]=iCell;   j_index[0]=jCell; /* cell (iCell,jCell) */
    /* First ring */
    i_index[1]=iCell-1; j_index[1]=jCell-1;
    i_index[2]=iCell;   j_index[2]=jCell-1;
    i_index[3]=iCell+1; j_index[3]=jCell-1;
    i_index[4]=iCell-1; j_index[4]=jCell;
    i_index[5]=iCell+1; j_index[5]=jCell;
    i_index[6]=iCell-1; j_index[6]=jCell+1;
    i_index[7]=iCell;   j_index[7]=jCell+1;
    i_index[8]=iCell+1; j_index[8]=jCell+1;
    break;

  default: // general expression
    i_index[0] = iCell;
    j_index[0] = jCell;
    for (int i=iCell-rings, Poz=1; i<=iCell+rings; ++i)
      for (int j=jCell-rings; j<=jCell+rings; ++j){
	if(!((i==iCell)&&(j==jCell)) ){
	  i_index[Poz] = i;
	  j_index[Poz] = j;
	  ++Poz;
	}
      }
  }//endswitch
}

// MakeReconstructionStencil(int,int,int,int,,vector<int>,vector<int>) function
// Set the stencil for the 2D kExact reconstruction used at curved boundaries
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       const int NorthCurvedBnd, const int SouthCurvedBnd,
			       const int EastCurvedBnd, const int WestCurvedBnd,
			       const int &ICl, const int &ICu, const int &JCl, const int &JCu,
			       int & StencilDimension, 
			       vector<int> & i_index, vector<int> & j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  i_index[0] = iCell;
  j_index[0] = jCell;
  int i,j, Imin, Imax, Jmin, Jmax, Poz;

  /* Determine Imin, Imax, Jmin, Jmax */
  Imin = iCell-rings; Imax = iCell+rings;
  Jmin = jCell-rings; Jmax = jCell+rings;

  if( NorthCurvedBnd==1 && Jmax > JCu ){ /* the north boundary is curved */
    if (jCell == JCu){
      Jmin -= 1;		/* add one extra layer in the j direction*/
      Jmax = JCu;
    } else {
      Jmax = JCu;		/* limit Jmax */
    }
  }

  if(SouthCurvedBnd==1 && Jmin < JCl){ 	/* the south boundary is curved */
    if (jCell == JCl){
      Jmax += 1;		/* add one extra layer in the j direction*/
      Jmin = JCl;
    } else {
      Jmin = JCl;		/* limit Jmin */
    }
  }

  if(EastCurvedBnd==1 && Imax > ICu){  /* the east boundary is curved */
    if (iCell == ICu){
      Imin -= 1;                /* add one extra layer in the i direction*/
      Imax = ICu;
    } else {
      Imax = ICu;		/* limit Imax */
    }
  }

  if(WestCurvedBnd==1 && Imin < ICl){  	/* the west boundary is curved */
    if (iCell == ICl){
      Imax += 1;                /* add one extra layer in the i direction*/
      Imin = ICl;
    } else {
      Imin = ICl;		/* limit Imin */
    }
  }

  /* Form stencil */
  for (i=Imin, Poz=1; i<=Imax; ++i)
    for (j=Jmin; j<=Jmax; ++j){
      if(!((i==iCell)&&(j==jCell)) ){
	i_index[Poz] = i;
	j_index[Poz] = j;
	++Poz;
      }//endif
    }// endfor

  StencilDimension = Poz;
}

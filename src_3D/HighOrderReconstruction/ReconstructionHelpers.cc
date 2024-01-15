/*!\file ReconstructionHelpers.cc
  \brief Source file providing implementation of subroutines prototyped in ReconstructionHelpers.h */

#include "ReconstructionHelpers.h"

/******************************************************************************************
Function for generating the geometric coefficients for a cartesian cell
******************************************************************************************/
double GeomCoeffCartesian(int p1, int p2, int p3, 
			  double deltaX, double deltaY, double deltaZ, 
			  double deltaXC, double deltaYC, double deltaZC){

  /* p1 -> the first power coefficient
     p2 -> the second power coefficient
     p3 -> the second power coefficient
     deltaX -> the grid size in the X direction
     deltaY -> the grid size in the Y direction
     deltaZ -> the grid size in the Z direction
     deltaXC -> the X distance between the center of the reconstructed cell and that of the cell used in the reconstruction
     deltaYC -> the Y distance between the center of the reconstructed cell and that of the cell used in the reconstruction
     deltaZC -> the Z distance between the center of the reconstructed cell and that of the cell used in the reconstruction

     Obs. To compute the coefficient of the reconstructed cell, deltaXC, deltaYC, and deltaZC must be ZERO. 
  */

  double val1, val2, val3;
  double coef_x1, coef_x2, coef_y1, coef_y2, coef_z1, coef_z2;;

  val1 = val2 = val3 = 0.0;
  coef_x1 = deltaX/2  + deltaXC;
  coef_x2 = -deltaX/2 + deltaXC;
  coef_y1 = deltaY/2  + deltaYC;
  coef_y2 = -deltaY/2 + deltaYC;
  coef_z1 = deltaZ/2  + deltaZC;
  coef_z2 = -deltaZ/2 + deltaZC;

  for (int m=1; m<=p1+1; m++){
    val1 += pow(coef_x1,p1+1-m)*pow(coef_x2,m-1);
  }
  for (int l=1; l<=p2+1; l++){
    val2 += pow(coef_y1,p2+1-l)*pow(coef_y2,l-1);
  }
  for (int n=1; n<=p3+1; n++){
    val3 += pow(coef_z1,p3+1-n)*pow(coef_z2,n-1);
  }

  return val1*val2*val3/((p1+1)*(p2+1)*(p3+1));
}

// MakeReconstructionStencil(int,int,int,vector<int>,vector<int>,vector<int>) function
// Set the stencil for the 3D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       vector<int> & i_index, vector<int> & j_index, vector<int> & k_index){

  // Obs. The first position (i_index[0],j_index[0],k_index[0]) corresponds to (iCell,jCell,kCell)
  i_index[0]=iCell;
  j_index[0]=jCell; 
  k_index[0]=kCell;     

  int Poz = 1;

  /* First layer */
  for (int k=kCell-1; k<=kCell+1; ++k){
    for (int j=jCell-1; j<=jCell+1; ++j){
      for (int i=iCell-1; i<=iCell+1; ++i){
        if( !(i==iCell && j==jCell && k==kCell) ){
          i_index[Poz] = i;
          j_index[Poz] = j;
          k_index[Poz] = k;
          ++Poz;      
        }
      } /*end for*/
    } /*end for*/
  } /*end for*/


  if (rings == 2){

    /* Second layer */
    for (int k=kCell-2; k<=kCell+2; ++k){
      for (int j=jCell-2; j<=jCell+2; ++j){
        for (int i=iCell-2; i<=iCell+2; ++i){
          // For k = -2 and k = 2 fill in all nine cells
          if( (k*k) > 1 ){
            i_index[Poz] = i;
            j_index[Poz] = j;
            k_index[Poz] = k;
            ++Poz;
          }
          // For k=-1, 0, or 1: fill in outer cells only,
          // since inner cells belong to first layer of stencil
          else if ( (i*i) > 1 || (j*j) > 1 ){
              i_index[Poz] = i;
              j_index[Poz] = j;
              k_index[Poz] = k;
              ++Poz;
          }
        } /*end for*/
      } /*end for*/
    } /*end for*/
  } /*end if*/

}

//MakeReconstructionStencil(int,int,int,int *,int *,int *) function
// Set the stencil for the 3D kExact reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       int *i_index, int *j_index, int *k_index){

  // Obs. The first position (i_index[0],j_index[0],k_index[0]) corresponds to (iCell,jCell,kCell)
  i_index[0]=iCell;
  j_index[0]=jCell; 
  k_index[0]=kCell;     

  int Poz = 1;

  /* First layer */
  for (int k=kCell-1; k<=kCell+1; ++k){
    for (int j=jCell-1; j<=jCell+1; ++j){
      for (int i=iCell-1; i<=iCell+1; ++i){
        if( !(i==iCell && j==jCell && k==kCell) ){
          i_index[Poz] = i;
          j_index[Poz] = j;
          k_index[Poz] = k;
          ++Poz;      
        }
      } /*end for*/
    } /*end for*/
  } /*end for*/


  if (rings == 2){

    /* Second layer */
    for (int k=kCell-2; k<=kCell+2; ++k){
      for (int j=jCell-2; j<=jCell+2; ++j){
        for (int i=iCell-2; i<=iCell+2; ++i){
          // For cells on level k = -2 and level k = 2 fill in all nine cells
          if( (k-kCell)*(k-kCell) > 1 ){
            i_index[Poz] = i;
            j_index[Poz] = j;
            k_index[Poz] = k;
            ++Poz;
          }
          // For cells on level k=-1, 0, or 1: fill in outer-ring (16) cells only,
          // since inner cells belong to the first layer of stencil
          else if ( (i-iCell)*(i-iCell) > 1 || (j-jCell)*(j-jCell) > 1 ){
              i_index[Poz] = i;
              j_index[Poz] = j;
              k_index[Poz] = k;
              ++Poz;
          }
        } /*end for*/
      } /*end for*/
    } /*end for*/
  } /*end if*/
}




// --> RR: ReconstructionHelpers.cc make stencil for curved boundaries
//// MakeReconstructionStencil(int,int,int,int,,vector<int>,vector<int>) function
//// Set the stencil for the 2D kExact reconstruction used at curved boundaries
//void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
//			       const int NorthCurvedBnd, const int SouthCurvedBnd,
//			       const int EastCurvedBnd, const int WestCurvedBnd,
//			       const int &ICl, const int &ICu, const int &JCl, const int &JCu,
//			       int & StencilDimension, 
//			       vector<int> & i_index, vector<int> & j_index){
//
//  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)
//
//  i_index[0] = iCell;
//  j_index[0] = jCell;
//  int i,j, Imin, Imax, Jmin, Jmax, Poz;
//
//  /* Determine Imin, Imax, Jmin, Jmax */
//  Imin = iCell-rings; Imax = iCell+rings;
//  Jmin = jCell-rings; Jmax = jCell+rings;
//
//  if( NorthCurvedBnd==1 && Jmax > JCu ){ /* the north boundary is curved */
//    if (jCell == JCu){
//      Jmin -= 1;		/* add one extra layer in the j direction*/
//      Jmax = JCu;
//    } else {
//      Jmax = JCu;		/* limit Jmax */
//    }
//  }
//
//  if(SouthCurvedBnd==1 && Jmin < JCl){ 	/* the south boundary is curved */
//    if (jCell == JCl){
//      Jmax += 1;		/* add one extra layer in the j direction*/
//      Jmin = JCl;
//    } else {
//      Jmin = JCl;		/* limit Jmin */
//    }
//  }
//
//  if(EastCurvedBnd==1 && Imax > ICu){  /* the east boundary is curved */
//    if (iCell == ICu){
//      Imin -= 1;                /* add one extra layer in the i direction*/
//      Imax = ICu;
//    } else {
//      Imax = ICu;		/* limit Imax */
//    }
//  }
//
//  if(WestCurvedBnd==1 && Imin < ICl){  	/* the west boundary is curved */
//    if (iCell == ICl){
//      Imax += 1;                /* add one extra layer in the i direction*/
//      Imin = ICl;
//    } else {
//      Imin = ICl;		/* limit Imin */
//    }
//  }
//
//  /* Form stencil */
//  for (i=Imin, Poz=1; i<=Imax; ++i)
//    for (j=Jmin; j<=Jmax; ++j){
//      if(!((i==iCell)&&(j==jCell)) ){
//	i_index[Poz] = i;
//	j_index[Poz] = j;
//	++Poz;
//      }//endif
//    }// endfor
//
//  StencilDimension = Poz;
//}

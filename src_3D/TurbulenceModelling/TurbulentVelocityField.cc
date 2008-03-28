
#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED
#include "TurbulentVelocityField.h"
#endif // _TURBULENT_VELOCITY_FIELD_INCLUDED


/***********************************************************************
 * Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction         *
 *                                 -- Reconstruct the velocity field.  *
 ***********************************************************************/
void Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction(const int &i,
								 const int &j,
								 const int &k,
								 Vector3D &dVdx,
								 Vector3D &dVdy,
								 Vector3D &dVdz) {
   int n, n2, n_pts, i_index[26], j_index[26], k_index[26];
   double DxDx_ave, DxDy_ave, DyDy_ave, DxDz_ave, DyDz_ave, DzDz_ave;
   double D;
   
   Vector3D  DU, DUDx_ave, DUDy_ave, DUDz_ave;
   Vector3D  D1, D2, D3;
   Vector3D  dX;


   if (i != ICl  &&  i != ICu &&
       j != JCl  &&  j != JCu &&
       k != KCl  &&  k != KCu) {

     n_pts = 26;
     // k plane
     i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
     i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
     i_index[2] = i+1; j_index[2] = j-1; k_index[2] = k;
     i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k;
     i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k;
     i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k;
     i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k;
     i_index[7] = i+1; j_index[7] = j+1; k_index[7] = k;
     //k-1 plane
     i_index[8] = i-1; j_index[8] = j-1; k_index[8] = k-1;
     i_index[9] = i  ; j_index[9] = j-1; k_index[9] = k-1;
     i_index[10] = i+1; j_index[10] = j-1; k_index[10] = k-1;
     i_index[11] = i-1; j_index[11] = j  ; k_index[11] = k-1;
     i_index[12] = i  ; j_index[12] = j  ; k_index[12] = k-1;
     i_index[13] = i+1; j_index[13] = j  ; k_index[13] = k-1;
     i_index[14] = i-1; j_index[14] = j+1; k_index[14] = k-1;
     i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k-1;
     i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k-1;
     //k+1 plane
     i_index[17] = i-1; j_index[17] = j-1; k_index[17] = k+1;
     i_index[18] = i  ; j_index[18] = j-1; k_index[18] = k+1;
     i_index[19] = i+1; j_index[19] = j-1; k_index[19] = k+1;
     i_index[20] = i-1; j_index[20] = j  ; k_index[20] = k+1;
     i_index[21] = i  ; j_index[21] = j  ; k_index[21] = k+1;
     i_index[22] = i+1; j_index[22] = j  ; k_index[22] = k+1;
     i_index[23] = i-1; j_index[23] = j+1; k_index[23] = k+1;
     i_index[24] = i  ; j_index[24] = j+1; k_index[24] = k+1;
     i_index[25] = i+1; j_index[25] = j+1; k_index[25] = k+1;

   } else {

     if (i == ICl) {

       if (j == JCl) {
	 // corner ICl, JCl, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   // k+1 plane
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k+1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k+1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k+1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k+1;

	   // corner ICl, JCl, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k-1;

	   // outer ICl, JCl
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane	 
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k-1;
	   //k+1 plane	 
	   i_index[7] = i  ; j_index[7] = j  ; k_index[7] = k+1;
	   i_index[8] = i+1; j_index[8] = j  ; k_index[8] = k+1;
	   i_index[9] = i  ; j_index[9] = j+1; k_index[9] = k+1;
	   i_index[10] = i+1; j_index[10] = j+1; k_index[10] = k+1;
	 }

       } else if (j == JCu) {
	 // corner ICl, JCu, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k+1 plane 
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k+1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k+1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k+1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k+1;
      
	   // corner ICl, JCu, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k-1;

	   // outer ICl, JCu
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i  ; j_index[7] = j-1; k_index[7] = k+1;
	   i_index[8] = i+1; j_index[8] = j-1; k_index[8] = k+1;
	   i_index[9] = i  ; j_index[9] = j  ; k_index[9] = k+1;
	   i_index[10] = i+1; j_index[10] = j  ; k_index[10] = k+1;
	 }

	 // inner ICl, != JCl, != JCu, != KCl, != KCu
       } else {
	 n_pts = 17;
	 // k plane
	 i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	 i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	 i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	 i_index[3] = i  ; j_index[3] = j+1; k_index[3] = k;
	 i_index[4] = i+1; j_index[4] = j+1; k_index[4] = k;
	 //k-1 plane
	 i_index[5] = i  ; j_index[5] = j-1; k_index[5] = k-1;
	 i_index[6] = i+1; j_index[6] = j-1; k_index[6] = k-1;
	 i_index[7] = i  ; j_index[7] = j  ; k_index[7] = k-1;
	 i_index[8] = i+1; j_index[8] = j  ; k_index[8] = k-1;
	 i_index[9] = i  ; j_index[9] = j+1; k_index[9] = k-1;
	 i_index[10] = i+1; j_index[10] = j+1; k_index[10] = k-1;
	 //k+1 plane
	 i_index[11] = i  ; j_index[11] = j-1; k_index[11] = k+1;
	 i_index[12] = i+1; j_index[12] = j-1; k_index[12] = k+1;
	 i_index[13] = i  ; j_index[13] = j  ; k_index[13] = k+1;
	 i_index[14] = i+1; j_index[14] = j  ; k_index[14] = k+1;
	 i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k+1;
	 i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k+1;
       }


     } else if (i == ICu) {

       if (j == JCl) {

	 // corner ICu, JCl, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k+1 plane
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k+1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k+1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k+1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k+1;

	   // corner ICu, JCl, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane      
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k-1;

	   // outer ICu, JCl
	 } else {
	   n_pts = 11;
	   // k plane 
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i-1; j_index[7] = j  ; k_index[7] = k+1;
	   i_index[8] = i  ; j_index[8] = j  ; k_index[8] = k+1;
	   i_index[9] = i-1; j_index[9] = j+1; k_index[9] = k+1;
	   i_index[10] = i  ; j_index[10] = j+1; k_index[10] = k+1;
	 }


       } else if (j == JCu) {
	 // corner ICu, JCu, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k+1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k+1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k+1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k+1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k+1; 

	   // corner ICu, JCu, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k-1;

	   // outer ICu, JCu
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i-1; j_index[7] = j-1; k_index[7] = k+1;
	   i_index[8] = i  ; j_index[8] = j-1; k_index[8] = k+1;      
	   i_index[9] = i-1; j_index[9] = j  ; k_index[9] = k+1;
	   i_index[10] = i  ; j_index[10] = j  ; k_index[10] = k+1;
	 }

	 // inner ICu, != JCl, != JCu, != KCl, != KCu
       } else {
	 n_pts = 17;
	 // k plane
	 i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	 i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	 i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	 i_index[3] = i-1; j_index[3] = j+1; k_index[3] = k;
	 i_index[4] = i  ; j_index[4] = j+1; k_index[4] = k;
	 //k-1 plane
	 i_index[5] = i-1; j_index[5] = j-1; k_index[5] = k-1;
	 i_index[6] = i  ; j_index[6] = j-1; k_index[6] = k-1;
	 i_index[7] = i-1; j_index[7] = j  ; k_index[7] = k-1;
	 i_index[8] = i  ; j_index[8] = j  ; k_index[8] = k-1;
	 i_index[9] = i-1; j_index[9] = j+1; k_index[9] = k-1;
	 i_index[10] = i  ; j_index[10] = j+1; k_index[10] = k-1;
	 //k+1 plane
	 i_index[11] = i-1; j_index[11] = j-1; k_index[11] = k+1;
	 i_index[12] = i  ; j_index[12] = j-1; k_index[12] = k+1;
	 i_index[13] = i-1; j_index[13] = j  ; k_index[13] = k+1;
	 i_index[14] = i  ; j_index[14] = j  ; k_index[14] = k+1;
	 i_index[15] = i-1; j_index[15] = j+1; k_index[15] = k+1;
	 i_index[16] = i  ; j_index[16] = j+1; k_index[16] = k+1;
       }

     }
    
   } /* end if */  
     

   if (n_pts > 0) {
      DUDx_ave.zero();
      DUDy_ave.zero();
      DUDz_ave.zero();
      D1.zero();
      D2.zero();
      D3.zero();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DxDz_ave = ZERO;
      DyDy_ave = ZERO;
      DyDz_ave = ZERO;
      DzDz_ave = ZERO;
      D = ZERO;
      
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
         dX = Position[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] - Position[i][j][k];
         DU = Velocity[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] - Velocity[i][j][k];
         
         DUDx_ave += DU*dX.x;
         DUDy_ave += DU*dX.y;
         DUDz_ave += DU*dX.z;
         DxDx_ave += dX.x*dX.x;
         DxDy_ave += dX.x*dX.y;
         DxDz_ave += dX.x*dX.z;
         DyDy_ave += dX.y*dX.y;
         DyDz_ave += dX.y*dX.z;
         DzDz_ave += dX.z*dX.z;
         
      } /* endfor */
      
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DUDz_ave = DUDz_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DxDz_ave = DxDz_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      DyDz_ave = DyDz_ave/double(n_pts);
      DzDz_ave = DzDz_ave/double(n_pts);
     
   
      // use cramer's rule for this simple system

      D = DxDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
          DxDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D1 = DUDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
           DUDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
           DUDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D2 =DxDx_ave*(DUDy_ave* DzDz_ave - DUDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DUDz_ave - DUDx_ave*DzDz_ave)+
          DxDz_ave*(DUDx_ave*DyDz_ave - DxDz_ave*DUDy_ave);

      D3 =DxDx_ave*(DyDy_ave* DUDz_ave - DyDz_ave*DUDy_ave) +
          DxDy_ave*(DUDx_ave*DyDz_ave - DxDy_ave*DUDz_ave)+
          DxDz_ave*(DxDy_ave*DUDy_ave - DUDx_ave*DyDy_ave);

      dVdx = D1/D;
      dVdy = D2/D;
      dVdz = D3/D;  

   } else {
      dVdx.zero();
      dVdy.zero();
      dVdz.zero();
   } /* endif */
    
}

/***********************************************************************
 * Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction         *
 *                                 -- Reconstruct the velocity field.  *
 ***********************************************************************/
void Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction(const int &i,
								 const int &j,
								 const int &k) {

   int n, n2, n_pts, i_index[26], j_index[26], k_index[26];
   double DxDx_ave, DxDy_ave, DyDy_ave, DxDz_ave, DyDz_ave, DzDz_ave;
   double D;
   
   Vector3D  DU, DUDx_ave, DUDy_ave, DUDz_ave;
   Vector3D  D1, D2, D3;
   Vector3D  dX;


   if (i != ICl  &&  i != ICu &&
       j != JCl  &&  j != JCu &&
       k != KCl  &&  k != KCu) {

     n_pts = 26;
     // k plane
     i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
     i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
     i_index[2] = i+1; j_index[2] = j-1; k_index[2] = k;
     i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k;
     i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k;
     i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k;
     i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k;
     i_index[7] = i+1; j_index[7] = j+1; k_index[7] = k;
     //k-1 plane
     i_index[8] = i-1; j_index[8] = j-1; k_index[8] = k-1;
     i_index[9] = i  ; j_index[9] = j-1; k_index[9] = k-1;
     i_index[10] = i+1; j_index[10] = j-1; k_index[10] = k-1;
     i_index[11] = i-1; j_index[11] = j  ; k_index[11] = k-1;
     i_index[12] = i  ; j_index[12] = j  ; k_index[12] = k-1;
     i_index[13] = i+1; j_index[13] = j  ; k_index[13] = k-1;
     i_index[14] = i-1; j_index[14] = j+1; k_index[14] = k-1;
     i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k-1;
     i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k-1;
     //k+1 plane
     i_index[17] = i-1; j_index[17] = j-1; k_index[17] = k+1;
     i_index[18] = i  ; j_index[18] = j-1; k_index[18] = k+1;
     i_index[19] = i+1; j_index[19] = j-1; k_index[19] = k+1;
     i_index[20] = i-1; j_index[20] = j  ; k_index[20] = k+1;
     i_index[21] = i  ; j_index[21] = j  ; k_index[21] = k+1;
     i_index[22] = i+1; j_index[22] = j  ; k_index[22] = k+1;
     i_index[23] = i-1; j_index[23] = j+1; k_index[23] = k+1;
     i_index[24] = i  ; j_index[24] = j+1; k_index[24] = k+1;
     i_index[25] = i+1; j_index[25] = j+1; k_index[25] = k+1;

   } else {

     if (i == ICl) {

       if (j == JCl) {
	 // corner ICl, JCl, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   // k+1 plane
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k+1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k+1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k+1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k+1;

	   // corner ICl, JCl, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k-1;

	   // outer ICl, JCl
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i+1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane	 
	   i_index[3] = i  ; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j+1; k_index[6] = k-1;
	   //k+1 plane	 
	   i_index[7] = i  ; j_index[7] = j  ; k_index[7] = k+1;
	   i_index[8] = i+1; j_index[8] = j  ; k_index[8] = k+1;
	   i_index[9] = i  ; j_index[9] = j+1; k_index[9] = k+1;
	   i_index[10] = i+1; j_index[10] = j+1; k_index[10] = k+1;
	 }

       } else if (j == JCu) {
	 // corner ICl, JCu, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k+1 plane 
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k+1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k+1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k+1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k+1;
      
	   // corner ICl, JCu, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k-1;

	   // outer ICl, JCu
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i  ; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i+1; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i  ; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i+1; j_index[6] = j  ; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i  ; j_index[7] = j-1; k_index[7] = k+1;
	   i_index[8] = i+1; j_index[8] = j-1; k_index[8] = k+1;
	   i_index[9] = i  ; j_index[9] = j  ; k_index[9] = k+1;
	   i_index[10] = i+1; j_index[10] = j  ; k_index[10] = k+1;
	 }

	 // inner ICl, != JCl, != JCu, != KCl, != KCu
       } else {
	 n_pts = 17;
	 // k plane
	 i_index[0] = i  ; j_index[0] = j-1; k_index[0] = k;
	 i_index[1] = i+1; j_index[1] = j-1; k_index[1] = k;
	 i_index[2] = i+1; j_index[2] = j  ; k_index[2] = k;
	 i_index[3] = i  ; j_index[3] = j+1; k_index[3] = k;
	 i_index[4] = i+1; j_index[4] = j+1; k_index[4] = k;
	 //k-1 plane
	 i_index[5] = i  ; j_index[5] = j-1; k_index[5] = k-1;
	 i_index[6] = i+1; j_index[6] = j-1; k_index[6] = k-1;
	 i_index[7] = i  ; j_index[7] = j  ; k_index[7] = k-1;
	 i_index[8] = i+1; j_index[8] = j  ; k_index[8] = k-1;
	 i_index[9] = i  ; j_index[9] = j+1; k_index[9] = k-1;
	 i_index[10] = i+1; j_index[10] = j+1; k_index[10] = k-1;
	 //k+1 plane
	 i_index[11] = i  ; j_index[11] = j-1; k_index[11] = k+1;
	 i_index[12] = i+1; j_index[12] = j-1; k_index[12] = k+1;
	 i_index[13] = i  ; j_index[13] = j  ; k_index[13] = k+1;
	 i_index[14] = i+1; j_index[14] = j  ; k_index[14] = k+1;
	 i_index[15] = i  ; j_index[15] = j+1; k_index[15] = k+1;
	 i_index[16] = i+1; j_index[16] = j+1; k_index[16] = k+1;
       }


     } else if (i == ICu) {

       if (j == JCl) {

	 // corner ICu, JCl, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k+1 plane
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k+1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k+1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k+1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k+1;

	   // corner ICu, JCl, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane      
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k-1;

	   // outer ICu, JCl
	 } else {
	   n_pts = 11;
	   // k plane 
	   i_index[0] = i-1; j_index[0] = j  ; k_index[0] = k;
	   i_index[1] = i-1; j_index[1] = j+1; k_index[1] = k;
	   i_index[2] = i  ; j_index[2] = j+1; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j  ; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j  ; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j+1; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j+1; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i-1; j_index[7] = j  ; k_index[7] = k+1;
	   i_index[8] = i  ; j_index[8] = j  ; k_index[8] = k+1;
	   i_index[9] = i-1; j_index[9] = j+1; k_index[9] = k+1;
	   i_index[10] = i  ; j_index[10] = j+1; k_index[10] = k+1;
	 }


       } else if (j == JCu) {
	 // corner ICu, JCu, KCl
	 if (k == KCl) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k+1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k+1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k+1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k+1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k+1; 

	   // corner ICu, JCu, KCu
	 } else if (k == KCu) {
	   n_pts = 7;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k-1;

	   // outer ICu, JCu
	 } else {
	   n_pts = 11;
	   // k plane
	   i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	   i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	   i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	   //k-1 plane
	   i_index[3] = i-1; j_index[3] = j-1; k_index[3] = k-1;
	   i_index[4] = i  ; j_index[4] = j-1; k_index[4] = k-1;
	   i_index[5] = i-1; j_index[5] = j  ; k_index[5] = k-1;
	   i_index[6] = i  ; j_index[6] = j  ; k_index[6] = k-1;
	   //k+1 plane
	   i_index[7] = i-1; j_index[7] = j-1; k_index[7] = k+1;
	   i_index[8] = i  ; j_index[8] = j-1; k_index[8] = k+1;      
	   i_index[9] = i-1; j_index[9] = j  ; k_index[9] = k+1;
	   i_index[10] = i  ; j_index[10] = j  ; k_index[10] = k+1;
	 }

	 // inner ICu, != JCl, != JCu, != KCl, != KCu
       } else {
	 n_pts = 17;
	 // k plane
	 i_index[0] = i-1; j_index[0] = j-1; k_index[0] = k;
	 i_index[1] = i  ; j_index[1] = j-1; k_index[1] = k;
	 i_index[2] = i-1; j_index[2] = j  ; k_index[2] = k;
	 i_index[3] = i-1; j_index[3] = j+1; k_index[3] = k;
	 i_index[4] = i  ; j_index[4] = j+1; k_index[4] = k;
	 //k-1 plane
	 i_index[5] = i-1; j_index[5] = j-1; k_index[5] = k-1;
	 i_index[6] = i  ; j_index[6] = j-1; k_index[6] = k-1;
	 i_index[7] = i-1; j_index[7] = j  ; k_index[7] = k-1;
	 i_index[8] = i  ; j_index[8] = j  ; k_index[8] = k-1;
	 i_index[9] = i-1; j_index[9] = j+1; k_index[9] = k-1;
	 i_index[10] = i  ; j_index[10] = j+1; k_index[10] = k-1;
	 //k+1 plane
	 i_index[11] = i-1; j_index[11] = j-1; k_index[11] = k+1;
	 i_index[12] = i  ; j_index[12] = j-1; k_index[12] = k+1;
	 i_index[13] = i-1; j_index[13] = j  ; k_index[13] = k+1;
	 i_index[14] = i  ; j_index[14] = j  ; k_index[14] = k+1;
	 i_index[15] = i-1; j_index[15] = j+1; k_index[15] = k+1;
	 i_index[16] = i  ; j_index[16] = j+1; k_index[16] = k+1;
       }

     }
    
   } /* end if */  
     

   if (n_pts > 0) {
      DUDx_ave.zero();
      DUDy_ave.zero();
      DUDz_ave.zero();
      D1.zero();
      D2.zero();
      D3.zero();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DxDz_ave = ZERO;
      DyDy_ave = ZERO;
      DyDz_ave = ZERO;
      DzDz_ave = ZERO;
      D = ZERO;
      
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
         dX = Position[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] - Position[i][j][k];
         DU = Velocity[ i_index[n2] ][ j_index[n2] ][ k_index[n2] ] - Velocity[i][j][k];
         
         DUDx_ave += DU*dX.x;
         DUDy_ave += DU*dX.y;
         DUDz_ave += DU*dX.z;
         DxDx_ave += dX.x*dX.x;
         DxDy_ave += dX.x*dX.y;
         DxDz_ave += dX.x*dX.z;
         DyDy_ave += dX.y*dX.y;
         DyDz_ave += dX.y*dX.z;
         DzDz_ave += dX.z*dX.z;
         
      } /* endfor */
      
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DUDz_ave = DUDz_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DxDz_ave = DxDz_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      DyDz_ave = DyDz_ave/double(n_pts);
      DzDz_ave = DzDz_ave/double(n_pts);
     
   
      // use cramer's rule for this simple system

      D = DxDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
          DxDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D1 = DUDx_ave*(DyDy_ave* DzDz_ave - DyDz_ave*DyDz_ave) +
           DUDy_ave*(DxDz_ave*DyDz_ave - DxDy_ave*DzDz_ave)+
           DUDz_ave*(DxDy_ave*DyDz_ave - DxDz_ave*DyDy_ave);
      
      D2 =DxDx_ave*(DUDy_ave* DzDz_ave - DUDz_ave*DyDz_ave) +
          DxDy_ave*(DxDz_ave*DUDz_ave - DUDx_ave*DzDz_ave)+
          DxDz_ave*(DUDx_ave*DyDz_ave - DxDz_ave*DUDy_ave);

      D3 =DxDx_ave*(DyDy_ave* DUDz_ave - DyDz_ave*DUDy_ave) +
          DxDy_ave*(DUDx_ave*DyDz_ave - DxDy_ave*DUDz_ave)+
          DxDz_ave*(DxDy_ave*DUDy_ave - DUDx_ave*DyDy_ave);

      dVdx[i][j][k] = D1/D;
      dVdy[i][j][k] = D2/D;
      dVdz[i][j][k] = D3/D;  

   } else {
      dVdx[i][j][k].zero();
      dVdy[i][j][k].zero();
      dVdz[i][j][k].zero();
   } /* endif */
    
}

/***********************************************************************
 * Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction         *
 *                                 -- Reconstruct the velocity field.  *
 ***********************************************************************/
void Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction(void) {
  for (int i = ICl; i <= ICu; ++i ) {
    for (int j  = JCl; j <= JCu; ++j ) {
      for (int  k  = KCl; k <= KCu; ++k ) {
	LeastSquares_Reconstruction(i, j, k);
      }
    }
  }   
      
}

/*********************************************************************
 * Turbulent_Velocity_Field_Block::Broadcast                         *
 *                     -- Broadcast turbulent velocity field block.  *
 *********************************************************************/
void Turbulent_Velocity_Field_Block::Broadcast(void) {

#ifdef _MPI_VERSION

  int ni, nj, nk, ng, block_allocated, buffer_size;
  double *buffer;

  int N_3Dvectors = 2; // Position and Velocity

  /* Broadcast the number of cells in each direction. */

  if (CFFC_Primary_MPI_Processor()) {
    ni = NCi;
    nj = NCj;
    nk = NCk;
    ng = Nghost;

    if (Allocated)
      block_allocated = 1;
    else
      block_allocated = 0;
     
  } /* endif */

  MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&nk, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&block_allocated, 1, MPI::INT, 0);

  /* On non-primary MPI processors, allocate (re-allocate) 
     memory for the block as necessary. */

  if (!CFFC_Primary_MPI_Processor()) {
    if (NCi != ni || NCj != nj ||  NCk != nk || Nghost != ng ) { 
      if (block_allocated) allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng); 
    }
  } 

  /* Broadcast the global block number and diagonally opposite corners. */
    
  MPI::COMM_WORLD.Bcast(&gblknum, 1, MPI::INT, 0);
  
  MPI::COMM_WORLD.Bcast(&(Node_INl_JNl_KNl.x),
                          1,
			  MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(Node_INl_JNl_KNl.y),
                          1,
			  MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(Node_INl_JNl_KNl.z),
                          1,
			  MPI::DOUBLE, 0);

  MPI::COMM_WORLD.Bcast(&(Node_INu_JNu_KNu.x),
                          1,
			  MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(Node_INu_JNu_KNu.y),
                          1,
			  MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&(Node_INu_JNu_KNu.z),
                          1,
			  MPI::DOUBLE, 0);
  
   
  /* Broadcast the velocities */

  if (block_allocated) {
    ni = (ICu+Nghost) - (ICl-Nghost) + 1;
    nj = (JCu+Nghost) - (JCl-Nghost) + 1;
    nk = (KCu+Nghost) - (KCl-Nghost) + 1;
    
    buffer_size = N_3Dvectors*(3*ni*nj*nk);
    buffer = new double[buffer_size]; 
    
    if (CFFC_Primary_MPI_Processor()) {
      int i_buffer = 0;
      for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	    for (int l=1;  l<=N_3Dvectors; ++l) {
	      for ( int n=1; n<=3; ++n) {
		if (l==1) {
		  buffer[i_buffer] = Velocity[i][j][k][n];
		} else if (l==2) {
		  buffer[i_buffer] = Position[i][j][k][n];
		}
		i_buffer ++;
	      }
	    }
	  }
	}
      } /* endfor */ 
    } /* endif */
	    
    MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

    if (!CFFC_Primary_MPI_Processor()) {
      int i_buffer = 0;
      for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
	    for (int l=1;  l<=N_3Dvectors; ++l) {
	      for ( int n=1; n<=3; ++n) {
		if (l==1) {
		  Velocity[i][j][k][n] = buffer[i_buffer];
		} else if (l==2) {
		  Position[i][j][k][n] = buffer[i_buffer];
		}
		i_buffer ++;
	      }
	    }
	  }	 
	}
      } /* endfor */
    } /* endif */
       
    delete []buffer; 
    buffer = NULL;

  } /* endif */
#endif

}


/*********************************************************************
 * Turbulent_Velocity_Field_Block::Broadcast                         *
 *                     -- Broadcast turbulent velocity field block.  *
 *********************************************************************/
#ifdef _MPI_VERSION
void Turbulent_Velocity_Field_Block::Broadcast(MPI::Intracomm &Communicator, 
					       const int Source_CPU){
    
  int Source_Rank = 0;
  int ni, nj, nk, ng, block_allocated, buffer_size;
  double *buffer;
    
  int N_3Dvectors = 2; // Position and Velocity
    
  /* Broadcast the number of cells in each direction. */
    
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    ni = NCi;
    nj = NCj;
    nk = NCk;
    ng = Nghost; 

    if (Allocated)
      block_allocated = 1;
    else
      block_allocated = 0;
  } /* endif */
    
  Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank); 
  Communicator.Bcast(&nk, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
  Communicator.Bcast(&block_allocated, 1, MPI::INT, Source_Rank);

    
  /* On non-source MPI processors, allocate (re-allocate) 
     memory for the block as necessary. */
    
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    if (NCi != ni || NCj != nj ||  NCk != nk || Nghost != ng) { 
      if (block_allocated) {
	allocate(ni-2*ng, nj-2*ng, nk-2*ng, ng); 
      }
    } /* endif */
  } /* endif */
    

  /* Broadcast the global block number and diagonally opposite corners. */
    
  Communicator.Bcast(&gblknum, 1, MPI::INT, Source_Rank);

  Communicator.Bcast(&(Node_INl_JNl_KNl.x),
		     1,
		     MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(Node_INl_JNl_KNl.y),
		     1,
		     MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(Node_INl_JNl_KNl.z),
		     1,
		     MPI::DOUBLE, Source_Rank);

  Communicator.Bcast(&(Node_INu_JNu_KNu.x),
		     1,
		     MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(Node_INu_JNu_KNu.y),
		     1,
		     MPI::DOUBLE, Source_Rank);
  Communicator.Bcast(&(Node_INu_JNu_KNu.z),
		     1,
		     MPI::DOUBLE, Source_Rank);
    
  /* Broadcast the velocities */
    
  if (block_allocated) {
    ni = (ICu+Nghost) - (ICl-Nghost) + 1;
    nj = (JCu+Nghost) - (JCl-Nghost) + 1;
    nk = (KCu+Nghost) - (KCl-Nghost) + 1;

    buffer_size = N_3Dvectors*(3*ni*nj*nk);
    buffer = new double[buffer_size];

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      int i_buffer = 0;
      for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
	    for (int l=1;  l<=N_3Dvectors; ++l) {
	      for ( int n=1; n<=3; ++n) {
		if (l==1) {
		  buffer[i_buffer] = Velocity[i][j][k][n];
		} else if (l==2) {
		  buffer[i_buffer] = Position[i][j][k][n];
		}
		i_buffer ++;
	      }
	    }
	  }	 
	}
      } /* endfor */
    } /* endif */
        
    Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      int i_buffer = 0;
      for (int k  = KCl-Nghost ; k <= KCu+Nghost ; ++k ) {
	for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
	  for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i) {
	    for (int l=1;  l<=N_3Dvectors; ++l) {
	      for ( int n=1; n<=3; ++n) {
		if (l==1) {
		  Velocity[i][j][k][n] = buffer[i_buffer];
		} else if (l==2) {
		  Position[i][j][k][n] = buffer[i_buffer];
		}                       
		i_buffer ++;
	      }
	    }
	  }	 
	}
      } /* endfor */
    } /* endif */

    delete []buffer; 
    buffer = NULL;
                
  } /* endif */ 

}
#endif 


/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Interpolate_Turbulent_Field    *
 *       -- Interpolate velocity field onto another grid.                    *
 *****************************************************************************/
void Turbulent_Velocity_Field_Multi_Block_List::
Interpolate_Turbulent_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			    Turbulent_Velocity_Field_Multi_Block_List &Interpolated_Velocity_Field) {

  int nnBlk, ii, jj, kk, n;
  double delta, dmin;
  double xmax, xmin, ymax, ymin, zmax, zmin;
  Vector3D dVdx, dVdy, dVdz, dX;

  for (int nBlk = 0; nBlk <= Initial_Mesh.NBlk-1; ++nBlk ) {
    for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
      for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
	for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {

	  Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k] = Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc;

	  //---------------------------
	  //   Cartesian box
	  //---------------------------

	  // find  nnBlk, ii, jj and kk to perform the interpolation/reconstruction
	  for (nnBlk = 0; nnBlk < NBlk; ++nnBlk) {
	    xmax = max(Vel_Blks[nnBlk].Node_INl_JNl_KNl.x, Vel_Blks[nnBlk].Node_INu_JNu_KNu.x);
	    xmin = min(Vel_Blks[nnBlk].Node_INl_JNl_KNl.x, Vel_Blks[nnBlk].Node_INu_JNu_KNu.x);

	    ymax = max(Vel_Blks[nnBlk].Node_INl_JNl_KNl.y, Vel_Blks[nnBlk].Node_INu_JNu_KNu.y);
	    ymin = min(Vel_Blks[nnBlk].Node_INl_JNl_KNl.y, Vel_Blks[nnBlk].Node_INu_JNu_KNu.y);

	    zmax = max(Vel_Blks[nnBlk].Node_INl_JNl_KNl.z, Vel_Blks[nnBlk].Node_INu_JNu_KNu.z);
	    zmin = min(Vel_Blks[nnBlk].Node_INl_JNl_KNl.z, Vel_Blks[nnBlk].Node_INu_JNu_KNu.z);

	    if ( (Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].x >= xmin  && 
		  Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].x <= xmax)  &&
		 (Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].y >= ymin  && 
		  Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].y <= ymax) &&
		 (Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].z >= zmin && 
		  Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].z <= zmax) ) {
	      break;
	    }

	  } /* end for*/   

	  
	  
	  // search in X-direction
	  dmin = 1E10;
	  for (n = Vel_Blks[nnBlk].ICl; n <= Vel_Blks[nnBlk].ICu; ++n) {
	    delta = fabs(Vel_Blks[nnBlk].Position[n][Vel_Blks[nnBlk].JCl][Vel_Blks[nnBlk].KCl].x -
			 Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].x);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      ii = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  // search in Y-direction
	  dmin = 1E10; 
	  for (n = Vel_Blks[nnBlk].JCl; n <= Vel_Blks[nnBlk].JCu; ++n) {
	    delta = fabs(Vel_Blks[nnBlk].Position[Vel_Blks[nnBlk].ICl][n][Vel_Blks[nnBlk].KCl].y -
			 Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].y);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      jj = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  // search in Z-direction
	  dmin = 1E10; 
	  for (n = Vel_Blks[nnBlk].KCl; n <= Vel_Blks[nnBlk].KCu; ++n) {
	    delta = fabs(Vel_Blks[nnBlk].Position[Vel_Blks[nnBlk].ICl][Vel_Blks[nnBlk].JCl][n].z - 
			 Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].z);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      kk = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  if (ii < Vel_Blks[nnBlk].ICl || ii > Vel_Blks[nnBlk].ICu || 
	      jj < Vel_Blks[nnBlk].JCl || jj > Vel_Blks[nnBlk].JCu || 
	      kk < Vel_Blks[nnBlk].KCl || kk > Vel_Blks[nnBlk].KCu) {
	    cout << "\n Index out of bound!!! -> ii = " << ii 
		 << "  jj = " << jj << "  kk = " << kk
		 << "\nICu = " << Vel_Blks[nnBlk].ICu 
		 << "  JCu = " << Vel_Blks[nnBlk].JCu 
		 << "  KCu = " << Vel_Blks[nnBlk].KCu; 
	  }


	  // use least squares to reconstruct the turbulent velocity field 
	  Vel_Blks[nnBlk].LeastSquares_Reconstruction(ii, jj, kk, 
						      dVdx, dVdy, dVdz);

	  dX = Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k] - 
	    Vel_Blks[nnBlk].Position[ii][jj][kk];

	  Interpolated_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k] = 
	    Vel_Blks[nnBlk].Velocity[ii][jj][kk] + dVdx*dX.x + dVdy*dX.y + dVdz*dX.z;

	} /* endfor */
      } /* endfor */
    } /* endfor */
  } /* endfor */

}



template <typename HEXA_BLOCK>
int Longitudinal_Correlation(const Octree_DataStructure &OcTree,
                             const AdaptiveBlock3D_ResourceList &Global_Soln_Block_List,
			     const AdaptiveBlock3D_List &LocalSolnBlockList,
			     HEXA_BLOCK *Solution_Block,
			     const Grid3D_Input_Parameters &IPs,
			     const double &u_ave, 
                             const double &sqr_u) {

  ofstream Out_Corr_Function;
  int error_flag = 0;
  int Nblks = Global_Soln_Block_List.Nused;
  double xref, yref, zref, uref, ucorr;
  double Yfuel_conditional = ZERO;
 
  int gblknum;
  double r, dr, max_r, fr, Lx, volume, local_vol, total_vol, R11, L11;
  Lx = IPs.Box_Width;
  dr = MILLION;

  //Conditional longitudinal correlation on fresh gas
  //  if (IPs.react_name != "NO_REACTIONS") {
  Yfuel_conditional = 0.95*0.05518;//IPs.Fresh_Fuel_Mass_Fraction;
    //  }

  for (int k = 0; k < LocalSolnBlockList.Nblk; ++k) {
    if (LocalSolnBlockList.Block[k].used == ADAPTIVEBLOCK3D_USED) {
      for (int ii = Solution_Block[k].ICl; ii <= Solution_Block[k].ICu; ++ii) {
        for (int jj = Solution_Block[k].JCl; jj <= Solution_Block[k].JCu; ++jj) {
           for (int kk = Solution_Block[k].KCl; kk <= Solution_Block[k].KCu; ++kk) {
	  if (Solution_Block[k].W[ii][jj][kk].spec[0].c >= Yfuel_conditional) {
	    dr = min(Solution_Block[k].Grid.Cell[ii][jj][kk].Xc.x - Solution_Block[k].Grid.Cell[ii-1][jj][kk].Xc.x, dr);
	  }
	}
      }
    }
   }
  }

  dr = CFFC_Minimum_MPI(dr);

  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.open("Correlation_Function.dat", ios::out);
    if(Out_Corr_Function.fail()){
      cerr<<"\n Error opening file: Correlation_Function.dat to write" << endl;
      exit(1);
    }
  }
  
  int *new_blocks_CPU, *CPUs_in_new_blocks;
  int my_rank, undefined_rank, new_blocks_BLK;
  CPUs_in_new_blocks = new int[Global_Soln_Block_List.Ncpu];
  new_blocks_CPU = new int[Global_Soln_Block_List.Ncpu];
  
#ifdef _MPI_VERSION
    MPI::Intracomm new_comm;
    MPI::Group     big_group = MPI::COMM_WORLD.Get_group();
    MPI::Group     new_group;
    undefined_rank = MPI::UNDEFINED;
#else
    undefined_rank = -1;
#endif
  
  HEXA_BLOCK  SolnBlk_duplicated;
  OctreeBlock *octree_block_duplicated_ptr;
  int SolnBlk_duplicated_info_level;

  Vector3D Vcorr, Xcorr;
  Vcorr.zero();  Xcorr.zero();
  bool correlated_flag, flag = false;
  int count = 0, count1 = 0;
  r = ZERO;
  L11 = ZERO;

/*   if (IPs.i_Grid == GRID_TURBULENT_PREMIXED_FLAME && */
/*       IPs.react_name == "NO_REACTIONS") { */
/*     max_r = HALF*(Lx-dr); */
/*   } else  { */
    max_r = Lx-dr;
/*   } */
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  while( r <= max_r) {
    total_vol = ZERO;
    R11 = ZERO;
    correlated_flag = false;
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    for (int iCPU = 0; iCPU <= OcTree.Ncpu-1; ++iCPU ) { // Loop over available processors.
      for (int iBLK = 0; iBLK <= OcTree.Nblk-1; ++iBLK ) { // Loop over available blocks.
        if (OcTree.Blocks[iCPU][iBLK] != NULL) {
          if (OcTree.Blocks[iCPU][iBLK]->block.used) { // Check if the solution block is used.
        
            octree_block_duplicated_ptr = OcTree.Blocks[iCPU][iBLK];
            new_blocks_CPU[0] = octree_block_duplicated_ptr->block.info.cpu;
            new_blocks_BLK = octree_block_duplicated_ptr->block.info.blknum;
            CPUs_in_new_blocks[0] = new_blocks_CPU[0];

	    if (Global_Soln_Block_List.Ncpu > 1) {
	      for (int iNEW = 1 ; iNEW < Global_Soln_Block_List.Ncpu; ++iNEW ) {
		if (iNEW != new_blocks_CPU[0]) {
		  new_blocks_CPU[iNEW] = iNEW;
		} else {
		  new_blocks_CPU[iNEW] = 0;
		}
		CPUs_in_new_blocks[iNEW] = new_blocks_CPU[iNEW];
	      }
	    }
  
#ifdef _MPI_VERSION
	    new_group = big_group.Incl(Global_Soln_Block_List.Ncpu, CPUs_in_new_blocks);
	    new_comm = MPI::COMM_WORLD.Create(new_group);
#endif
/* 	    if (LocalSolnBlockList.ThisCPU == new_blocks_CPU[0]) { */
/* 	      Copy_Solution_Block(SolnBlk_duplicated, Solution_Block[new_blocks_BLK]); */
/* 	      SolnBlk_duplicated_info_level = LocalSolnBlockList.Block[new_blocks_BLK].info.level; */
/* 	    } */
/* #ifdef _MPI_VERSION */
/* 	    if (my_rank != undefined_rank) { */
/* 	      Broadcast_Solution_Block(SolnBlk_duplicated, new_comm, new_blocks_CPU[0]); */
/* 	      new_comm.Bcast(&SolnBlk_duplicated_info_level, 1, MPI::INT, 0); */
/* 	    } */
/* #endif */

#ifdef _MPI_VERSION
	    if (new_comm != MPI::COMM_NULL) new_comm.Free();
	    new_group.Free();
#endif
	    for (int i_ref = SolnBlk_duplicated.ICl; i_ref <= SolnBlk_duplicated.ICu; ++i_ref) {
	      for (int j_ref = SolnBlk_duplicated.JCl; j_ref <= SolnBlk_duplicated.JCu; ++j_ref) {
  	         for (int k_ref = SolnBlk_duplicated.KCl; k_ref <= SolnBlk_duplicated.KCu; ++k_ref) {
/*                 if (SolnBlk_duplicated.W[i_ref][j_ref][k_ref].spec[0].c >= Yfuel_conditional) { */
/* 		  xref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.x; */
/* 		  yref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.y; */
/* 		  zref = SolnBlk_duplicated.Grid.Cell[i_ref][j_ref][k_ref].Xc.z; */
/* 		  uref = SolnBlk_duplicated.W[i_ref][j_ref][k_ref].v.x; */
/* 		  flag = false; */
/* 		} else { */
/*                   continue; */
/* 		} */
	       
		for (int q = 0; q < LocalSolnBlockList.Nblk; ++q) {
		  if (LocalSolnBlockList.Block[q].used == ADAPTIVEBLOCK3D_USED) {
		    for (int i = Solution_Block[q].ICl; i <= Solution_Block[q].ICu; ++i) {
		      for (int j = Solution_Block[q].JCl; j <= Solution_Block[q].JCu; ++j) {
		         for (int k = Solution_Block[q].KCl; k <= Solution_Block[q].KCu; ++k) {

			if ((Solution_Block[q].W[i][j][k].spec[0].c >= Yfuel_conditional) &&
                            (Solution_Block[q].Grid.xfaceW(i,j,k).x < (xref + r) &&
			     Solution_Block[q].Grid.xfaceE(i,j,k).x > (xref + r))) {

			  // Finer reference block than local block
			  if (SolnBlk_duplicated_info_level >= LocalSolnBlockList.Block[q].info.level  &&
			      (Solution_Block[q].Grid.xfaceS(i,j,k).y < yref  &&  Solution_Block[q].Grid.xfaceN(i,j,k).y > yref)) {

			    if (Solution_Block[q].Grid.Cell[i][j][k].Xc.y == yref &&
				Solution_Block[q].Grid.Cell[i][j][k].Xc.x == (xref + r)) {
			      Vcorr.x = Solution_Block[q].W[i][j][k].v.x;
			    } else {
			      Vector2D dX;
			      dX.x = (xref + r) - Solution_Block[q].Grid.Cell[i][j][k].Xc.x;
			      dX.y = yref - Solution_Block[q].Grid.Cell[i][j][k].Xc.y;
			      Vcorr.x = Solution_Block[q].W[i][j][k].v.x + Solution_Block[q].phi[i][j][k].v.x
				* (Solution_Block[q].dWdx[i][j][k].v.x * dX.x + Solution_Block[q].dWdy[i][j][k].v.x * dX.y);
			    }
			    // correlated u
			    ucorr = Vcorr.x;
			    local_vol = Solution_Block[q].Grid.volume(i,j,k);
			    count1 ++;
			    flag = true;

			  // Coarser reference block than local block
			  } else if (SolnBlk_duplicated_info_level < LocalSolnBlockList.Block[q].info.level  &&
				     (Solution_Block[q].Grid.Cell[i][j][k].Xc.y < yref  &&  Solution_Block[q].Grid.Cell[i][j+1][k].Xc.y > yref)) {
			    Xcorr.x = xref + r;
			    Xcorr.y = yref;
                            Trilinear_Interpolation(Solution_Block[q].Grid.Cell[i-1][j][k].Xc, Solution_Block[q].W[i-1][j][k],
                                                    Solution_Block[q].Grid.Cell[i][j][k].Xc, Solution_Block[q].W[i][j][k],
                                                    Solution_Block[q].Grid.Cell[i][j-1][k].Xc, Solution_Block[q].W[i][j-1][k],
                                                    Solution_Block[q].Grid.Cell[i-1][j-1][k].Xc, Solution_Block[q].W[i-1][j-1][k],
                                                    Solution_Block[q].Grid.Cell[i-1][j][k-1].Xc, Solution_Block[q].W[i-1][j][k-1],
                                                    Solution_Block[q].Grid.Cell[i][j][k-1].Xc, Solution_Block[q].W[i][j][k-1],
                                                    Solution_Block[q].Grid.Cell[i][j-1][k-1].Xc, Solution_Block[q].W[i][j-1][k-1],
                                                    Solution_Block[q].Grid.Cell[i-1][j-1][k-1].Xc, Solution_Block[q].W[i-1][j-1][k-1],
                                                    Solution_Block[q].Grid.Node[i][j][k].X);
/* 			    Bilinear_Interpolation(Solution_Block[q].W[i][j][k].v, Solution_Block[q].Grid.Cell[i][j][k].Xc, */
/* 						   Solution_Block[q].W[i][j+1][k].v, Solution_Block[q].Grid.Cell[i][j+1][k].Xc, */
/* 						   Solution_Block[q].W[i+1][j+1][k].v, Solution_Block[q].Grid.Cell[i+1][j+1][k].Xc, */
/* 						   Solution_Block[q].W[i+1][j][k].v, Solution_Block[q].Grid.Cell[i+1][j][k].Xc, */
/* 						   Xcorr, Vcorr); */
			    // correlated u
			    ucorr = Vcorr.x;
			    local_vol = Solution_Block[q].Grid.volume(i,j,k);
			    count1 ++;
			    flag = true;
			  } // end if
                                                                                          
// 			  if (Soln_ptr[q].Grid.Cell[i][j].Xc.y == yref &&
// 			      Soln_ptr[q].Grid.Cell[i][j].Xc.x == (xref + r)) {
// 			    Vcorr.x = Soln_ptr[q].W[i][j].v.x;
// 			  } else {
// 			    Xcorr.x = xref + r;
// 			    Xcorr.y = yref;
// 			    // Bilinear_Interpolation(Soln_ptr[q].WnNW(i,j).v, Soln_ptr[q].Grid.nodeNW(i,j).X,
// 			    // 						 Soln_ptr[q].WnNE(i,j).v, Soln_ptr[q].Grid.nodeNE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSE(i,j).v, Soln_ptr[q].Grid.nodeSE(i,j).X,
// 			    // 						 Soln_ptr[q].WnSW(i,j).v, Soln_ptr[q].Grid.nodeSW(i,j).X,
// 			    // 						 Xcorr, Vcorr);
// 			    Bilinear_Interpolation(Soln_ptr[q].W[i-1][j].v, Soln_ptr[q].Grid.Cell[i-1][j].Xc,
// 						   Soln_ptr[q].W[i][j+1].v, Soln_ptr[q].Grid.Cell[i][j+1].Xc,
// 						   Soln_ptr[q].W[i+1][j].v, Soln_ptr[q].Grid.Cell[i+1][j].Xc,
// 						   Soln_ptr[q].W[i][j-1].v, Soln_ptr[q].Grid.Cell[i][j-1].Xc,
// 						   Xcorr, Vcorr);
// 			  }
		  
			} // end if
		      }
		    }
		  } //end if
		  if (flag) break;
		} // end for
        
		if (flag) {
		  volume = SolnBlk_duplicated.Grid.volume(i_ref,j_ref,k_ref) + local_vol;
		  total_vol += volume;
		  R11 += (ucorr - u_ave)*(uref - u_ave)*volume;
		  correlated_flag = true;
		  count++;
		}
           
		} // end for
		 } // end for
	      }
	    }
	  }
	}
      }
    }
    if (CFFC_OR_MPI(correlated_flag)) {
      total_vol = CFFC_Summation_MPI(total_vol);
      R11 = CFFC_Summation_MPI(R11);
       
      R11 = R11/total_vol;     // Area-averaged two-point correlation
      fr = R11/sqr_u;        // Longitudinal autocorrelation function
      L11 += fr*dr;          // Integrate the above funtion to obtain L11
      if (CFFC_Primary_MPI_Processor()) {
        Out_Corr_Function << r << "  " << fr << endl;
      }
    }
       
    //cout << "\n->" << r;
    r += dr;
  } // end while
  
  count = CFFC_Summation_MPI(count);
  count1 = CFFC_Summation_MPI(count1);
    
  if (CFFC_Primary_MPI_Processor()) {
    Out_Corr_Function.close();
    cout << "\n\n *** L11 = " << L11 << " ***" << endl;
  }

  delete[] CPUs_in_new_blocks;   CPUs_in_new_blocks = NULL;
  delete[] new_blocks_CPU;       new_blocks_CPU = NULL;
  
  if (CPUs_in_new_blocks != NULL || new_blocks_CPU != NULL) error_flag = 1;
       
  return (error_flag);

}

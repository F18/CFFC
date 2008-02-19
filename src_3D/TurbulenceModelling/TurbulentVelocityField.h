#ifndef _TURBULENT_VELOCITY_FIELD_INCLUDED 
#define _TURBULENT_VELOCITY_FIELD_INCLUDED

/* Include required C++ header files. */

#include <cmath> 
#include <cassert>
#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function
#include <limits>
#include <complex>

/* Include required CFFC header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif // _VECTOR3D_INCLUDED

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif // _TENSOR3D_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "../Grid/Grid3DHexaMultiBlock.h"
#endif // _GRID3D_HEXA_MULTIBLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif // _FFTW_INCLUDED

// Constants
const complex<double>  I(0.0, 1.0);      // sqrt(-1.0)

#define TURBULENT_VELOCITY_FIELD_DATA_USED        1

#define TURBULENT_VELOCITY_FIELD_DATA_NOT_USED    0

/*!
 * Class: Turbulent_Velocity_Field_Block
 *
 * \brief Class used to store turbulent velocity field for a single block.
 *
 */
class Turbulent_Velocity_Field_Block {
  public:
    int             NCi,ICl,ICu; // i-direction turbulent velocity field cell counters
    int             NCj,JCl,JCu; // j-direction turbulent velocity field cell counters
    int             NCk,KCl,KCu; // k-direction turbulent velocity field cell counters
    int                  Nghost; // number of ghost cells
    Vector3D        ***Velocity; // array of turbulent velocity field vectors
    Vector3D        ***Position; // array of turbulent velocity field position vectors
    Vector3D   Node_INl_JNl_KNl, // diagonally opposite corners of the block 
               Node_INu_JNu_KNu;


    int Allocated; // Indicates whether or not the turbulent velocity field data has been allocated.
    
    /* Constructors. */
    Turbulent_Velocity_Field_Block(void) {
       NCi = 0; ICl = 0; ICu = 0; 
       NCj = 0; JCl = 0; JCu = 0;
       NCk = 0; KCl = 0; KCu = 0;
       Nghost = 0;
       Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
       Velocity = NULL;
       Position = NULL;
    }

    Turbulent_Velocity_Field_Block(const int Ni, 
                                   const int Nj, 
                                   const int Nk,
                                   const int Ng) {
      allocate(Ni, Nj, Nk, Ng);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Block(void) {
       deallocate();
    }

    /* Allocate memory for velocity field data. */
    void allocate(const int Ni, 
                  const int Nj, 
                  const int Nk,
                  const int Ng);
    
    /* Deallocate memory for velocity field data. */
    void deallocate(void);

    /* Reconstruct velocity field data. */
    void LeastSquares_Reconstruction(const int i,
				     const int j,
				     const int k,
				     Vector3D &dVdx,
				     Vector3D &dVdy,
				     Vector3D &dVdz);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, 
                                 const Turbulent_Velocity_Field_Block &V);

    friend istream &operator >> (istream &in_file, 
                                 Turbulent_Velocity_Field_Block &V);
  
    /* Other useful member functions. */

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Block(const Turbulent_Velocity_Field_Block &V);
    Turbulent_Velocity_Field_Block &operator =(const Turbulent_Velocity_Field_Block &V);
};

/*************************************************************************
 * Turbulent_Velocity_Field_Block::allocate -- Allocate memory.          *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::allocate(const int Ni, 
                                                     const int Nj, 
                                                     const int Nk,
                                                     const int Ng) {
   assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && Ng >=1 && !Allocated);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1;
   NCk = Nk+2*Ng; KCl = Ng; KCu = Nk+Ng-1;
   Nghost = Ng;
   Allocated = TURBULENT_VELOCITY_FIELD_DATA_USED;

   Velocity = new Vector3D**[NCi];
   Position = new Vector3D**[NCi];
   for (int i = 0; i <= NCi-1; ++i ){
      Velocity[i] = new Vector3D*[NCj];
      Position[i] = new Vector3D*[NCj];
      for (int j = 0; j <= NCj-1; ++j ){
         Velocity[i][j] = new Vector3D[NCk];
	 Position[i][j] = new Vector3D[NCk];
      } /* endfor */
   } /* endfor */
}

/*************************************************************************
 * Turbulent_Velocity_Field_Block::deallocate -- Deallocate memory.      *
 *************************************************************************/
inline void Turbulent_Velocity_Field_Block::deallocate(void) {
   if (Allocated) {
      assert(NCi >= 1 && NCj >= 1 && NCk >= 1);
      for (int i = 0; i <= NCi-1 ; ++i ) {
         for ( int j = 0 ; j <= NCj-1 ; ++j) {
            delete []Velocity[i][j]; Velocity[i][j] = NULL;
	    delete []Position[i][j]; Position[i][j] = NULL;
         } /* endfor */
         delete []Velocity[i]; Velocity[i] = NULL;
	 delete []Position[i]; Position[i] = NULL;
      }/*endfor*/
      delete []Velocity; Velocity = NULL;
      delete []Position; Position = NULL;
  
      NCi = 0; ICl = 0; ICu = 0; 
      NCj = 0; JCl = 0; JCu = 0; 
      NCk = 0; KCl = 0; KCu = 0; 
      Nghost = 0;
      Allocated = TURBULENT_VELOCITY_FIELD_DATA_NOT_USED;
   } /* endif */
}

/***********************************************************************
 * Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction         *
 *                                 -- Reconstruct the velocity field.  *
 ***********************************************************************/
inline void Turbulent_Velocity_Field_Block::LeastSquares_Reconstruction(const int i,
									const int j,
									const int k,
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



/*!
 * Class: Turbulent_Velocity_Field_Multi_Block
 *
 * \brief Class used to store turbulent velocity field for a
 *        1D array of turbulent velocity field solution blocks.
 *
 */
class Turbulent_Velocity_Field_Multi_Block_List {
  public:
    Turbulent_Velocity_Field_Block  *Vel_Blks; // one dimensional array of velocity block.
    int                                  NBlk;
    int                             NBlk_Idir, 
                                    NBlk_Jdir, 
                                    NBlk_Kdir; // Number of blocks in i, j and k directions.
    int                             Allocated; // Indicates if the velocity blocks have been allocated or not.
 
    /* Creation constructors. */
    Turbulent_Velocity_Field_Multi_Block_List(void) : 
       NBlk(0), NBlk_Idir(0), NBlk_Jdir(0), NBlk_Kdir(0), Vel_Blks(NULL), Allocated(0) { }

    Turbulent_Velocity_Field_Multi_Block_List(const int N) {
       Allocate(N);
    }

    Turbulent_Velocity_Field_Multi_Block_List(const int Ni, 
                                              const int Nj, 
                                              const int Nk) {
       Allocate(Ni, Nj, Nk);
    }

    /* Destructor. */
    ~Turbulent_Velocity_Field_Multi_Block_List(void) {
        Deallocate();
    }

    /* Other member functions  */

    void Allocate(const int Ni, const int Nj, const int Nk);

    void Allocate(const int N);

    void Deallocate(void);

    void Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                const Grid3D_Input_Parameters &Input);

    void Interpolate_Turbulent_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
				     Turbulent_Velocity_Field_Multi_Block_List &Interpolated_Velocity_Field);
    

  private:
    //copy and assignment are not permitted
    Turbulent_Velocity_Field_Multi_Block_List(const Turbulent_Velocity_Field_Multi_Block_List &V);
    Turbulent_Velocity_Field_Multi_Block_List &operator = (const Turbulent_Velocity_Field_Multi_Block_List &V);
};

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int Ni, 
                                                                const int Nj, 
                                                                const int Nk) {
   if (Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated) {
      NBlk_Idir = Ni; 
      NBlk_Jdir = Nj; 
      NBlk_Kdir = Nk; 
      NBlk = Ni*Nj*Nk;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Allocate memory.   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Allocate(const int N) {
   if (N >= 1 && !Allocated) {
      NBlk_Idir = N; 
      NBlk_Jdir = 1; 
      NBlk_Kdir = 1; 
      NBlk = N;
      Vel_Blks = new Turbulent_Velocity_Field_Block[NBlk];
      Allocated = 1;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Allocate -- Deallocate memory. *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Deallocate(void) {
   if (NBlk >= 1 && Allocated) {
       delete []Vel_Blks;
       Vel_Blks = NULL;
       NBlk_Idir = 0; 
       NBlk_Jdir = 0; 
       NBlk_Kdir = 0;
       NBlk = 0;
       Allocated = 0;
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Create --                      *
 *       Create memory containers for velocity field data.                   *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::Create(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                              const Grid3D_Input_Parameters &Input) {
   if (Initial_Mesh.NBlk >= 1 && Initial_Mesh.Allocated) {
      Allocate(Initial_Mesh.NBlk_Idir, Initial_Mesh.NBlk_Jdir, Initial_Mesh.NBlk_Kdir);
      for (int nBlk = 0; nBlk <= NBlk-1; ++nBlk ) {
	if (Initial_Mesh.Grid_Blks[nBlk].Allocated) Vel_Blks[nBlk].allocate(Initial_Mesh.Grid_Blks[nBlk].NCi-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCj-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].NCk-
                                                                            2*Initial_Mesh.Grid_Blks[nBlk].Nghost,
                                                                            Initial_Mesh.Grid_Blks[nBlk].Nghost);
      } /* endfor */
   } /* endif */
}

/*****************************************************************************
 * Turbulent_Velocity_Field_Multi_Block_List::Interpolate_Turbulent_Field    *
 *       -- Interpolate velocity field onto another grid.                    *
 *****************************************************************************/
inline void Turbulent_Velocity_Field_Multi_Block_List::
Interpolate_Turbulent_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			    Turbulent_Velocity_Field_Multi_Block_List &Interpolated_Velocity_Field) {

  int nnBlk, ii, jj, kk, n;
  double delta, dmin;
  double xmax, xmin, ymax, ymin, zmax, zmin;
  Vector3D dVdx, dVdy, dVdz, dX;

  int counter(0);  

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

	      counter++;
	      break;
	    }

	  } /* end for*/   

	  
	  
	  // search in X-direction
	  dmin = 1E9;
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
	  dmin = 1E9; 
	  for (n = Vel_Blks[nnBlk].JCl; n <= Vel_Blks[nnBlk].JCu; ++n) {
	    delta = fabs(Vel_Blks[nnBlk].Position[Vel_Blks[nnBlk].ICl][n][Vel_Blks[nnBlk].KCl].y -
			 Interpolated_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k].y);
	    if ( delta < dmin )  { 
	      dmin = delta;
	      jj = n;
	    }
	    if (dmin == 0.0) break;
	  } /* end for */


	  // search Z-direction
	  dmin = 1E9; 
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

  //cout << "\nInterpolation of turbulent field complete";
  //cout << "\n counter = " << counter; 

}



/*!
 * Class: RandomFieldRogallo
 *
 * \brief Class defined to generate random fluctuations using
 * Rogallo's procedure.
 *
 */
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
class RandomFieldRogallo {
  private:
    int spectrum_flag;  //!< Turbulence kinetic energy spectrum flag.

  public:

    //! Default constructor.
    RandomFieldRogallo() : spectrum_flag(SPECTRUM_VON_KARMAN_PAO) { }

    //! Another constructor.
    RandomFieldRogallo(const int &SPECTRUM) : spectrum_flag(SPECTRUM) { }  


    double k_1(const int &n1, const double &L1) const { return double(n1); }//2.0*PI*n1/L1; }//double(n1); }

    double k_2(const int &n2, const double &L2) const { return double(n2); }//2.0*PI*n2/L2; }//double(n2); }

    double k_3(const int &n3, const double &L3) const { return double(n3); }//2.0*PI*n3/L3; }//double(n3); }

    double random_double() const { return  drand48(); }

    // alpha and beta as defined by Rogallo, 1981
    complex<double> alpha(const double &abs_wave_num, 
			  const double &theta1, 
			  const double &phi) const;

    complex<double> beta(const double &abs_wave_num, 
			 const double &theta2, 
			 const double &phi) const;

    double Energy_Spectrum_Value(const double &abs_wave_num) const;

    int Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                                     Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
			                             const Grid3D_Input_Parameters &IPs);

    void Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                        const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
					const Grid3D_Input_Parameters &IPs) const;
       
    void Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                       Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
				       Grid3D_Input_Parameters &IPs);

};

//-----------------------------------------------------//
//         Members of RandomFieldRogallo class         //
//-----------------------------------------------------//

// alpha
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha(const double &abs_wave_num, const double &theta1, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta1) * cos(phi);
}

// beta
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta(const double &abs_wave_num, const double &theta2, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta2) * sin(phi);
}

// Prescribed energy spectrum
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Energy_Spectrum_Value(const double &abs_wave_num) const {

  double k, kp, kd, eps, u, Lp, EE;
  double C, A, alpha, a_s;
  int s;
 
  k = abs_wave_num;
  
  switch (spectrum_flag) {
    /*****  Lee and Reynolds  *****/
  case SPECTRUM_LEE_REYNOLDS :
    C = 20.0;  kp = 80.0;   
    EE = (k <= kp)  ?  C*k*k : C*kp*kp*pow(k/kp, -5.0/3.0);
    break;

    /*****  Laval and Nazarenko paper  *****/
  case SPECTRUM_LAVAL_NAZARENKO :
    //double EE, kp, C;
    C = 1.0;  kp = 4.0;
    EE = C*k*exp(-pow(k/kp, 2.0));
    break;

    /*****   von Karman-Pao   *****/
  case SPECTRUM_VON_KARMAN_PAO :
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp =1.0;   kd = 300.0 /*362.0*/;
    u = 0.8;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****  Haworth and Poinsot paper  *****/
  case SPECTRUM_HAWORTH_POINSOT :
    //double EE, kp, 
    u = 7.1;  Lp = TWO*PI/6.0;  // u=14.1
    kp = TWO*PI/Lp;
    EE = (32.0/3.0) * sqrt(2.0/PI)* (u*u/kp) * pow(k/kp, 4.0) * exp(-2.0*(k/kp)*(k/kp));
    break;

    /*****  Chasnov paper 1996  *****/
  case SPECTRUM_CHASNOV :  
    //double a_s, EE, kp = 4.0 /*8.0 20.0 4.0*/, u = 0.095; /*0.1 28.3 21.21  0.001*/  
    //int s = 3;
    s = 3;
    kp = 4.0;
    u = 0.009; 
    // u = 0.001  ->  Re_lambda = 10 
    // u = 0.009  ->  Re_lambda = 98
    // u = 0.095  ->  Re_lambda = 950
   
    a_s = pow(2.0*double(s)+1.0, double(s)+1.0)/(factorial(s)*pow(2.0, double(s)));
    EE = (HALF*a_s*u*u/kp)*pow(k/kp, 2.0*double(s)+1.0);
    EE = EE*exp(-(double(s)+HALF)*pow(k/kp, 2.0));
    break;

    /*****   Bell & Day report   *****/
  case SPECTRUM_BELL_DAY :
    //  kd = 1/(2*dx)
    kp = 3.0;   kd = ONE/0.576E-3;
    EE = pow(k/kp, 4.0) * exp(-9.0*pow(k/kd, 4.0/3.0)/4.0);
    EE /= pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****   von Karman-Pao   *****/
  default :
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp = 3.0;   kd = 1000.0 /*362.0*/;
    u = 0.107;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;    

  }  // end switch

  return (k == 0.0)  ?  0.0 : EE;
 
}
 
// Create_Homogeneous_Turbulence_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Create_Homogeneous_Turbulence_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
                                             Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                                             const Grid3D_Input_Parameters &IPs) {

  double L1, L2, L3;
  L1 = IPs.Box_Length;
  L2 = IPs.Box_Width;
  L3 = IPs.Box_Height;
  
  int Nx, Ny, Nz;
  Nx = IPs.NCells_Turbulence_Idir*Initial_Mesh.NBlk_Idir;
  Ny = IPs.NCells_Turbulence_Jdir*Initial_Mesh.NBlk_Jdir;
  Nz = IPs.NCells_Turbulence_Kdir*Initial_Mesh.NBlk_Kdir;

  
  double scaling_factor = 1.0/double(Nx*Ny*Nz);  // Scaling factor for the complex to real transform

  double        *u, *v, *w;        // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv, *ww;     // Arrays to store the velocity fluctuations in spectral space
  fftw_plan      physical;

  int index;
  int nz = Nz/2+1;

  // Allocation of arrays used in the transforms
  u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  w = (double *) malloc(Nx*Ny*Nz * sizeof(double));

  uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  
  int seed = 1; 
  //int seed = time(NULL);   // assigns the current time to the seed
  srand48(seed);             // changes the seed for drand48()
  int iconj, jconj;          // Position of the conjugate complex for the i and j indeces
  double k1, k2, k3;         // Wave numbers
  
  double theta1, theta2, phi;
  complex<double> aa, bb;
  double deno;

  if (CFFC_Primary_MPI_Processor()) {
    cout << "\n\n ==========================================================================\n"; 
    cout << " Generating Homogeneous Isotropic Turbulent Velocity Field"<<endl;
    cout << " ==========================================================================";
  }

  for (int i=0; i<Nx; ++i) {
    iconj = (i==0  ?  0 : Nx-i);
    for (int j=0; j<Ny; ++j) {
      jconj = (j==0  ?  0 : Ny-j);
      for (int l=0; l<Nz/2+1; ++l) {
	//Components of the wave number vector
	if( i<=Nx/2) {
	  k1 = k_1(i, L1);
	} else {
	  k1 = k_1(i-Nx, L1);
	}

	if( j<=Ny/2) {
	  k2 = k_2(j, L2);
	} else {
	  k2 = k_2(j-Ny, L2);
	}

	k3 = k_3(l, L3);

	// Wave number magnitude
	double abs_k = sqrt(k1*k1 + k2*k2 + k3*k3);

	//Rogallo's function  
	theta1 = 2.0*PI*random_double();  // Random number (0, 2*PI)
	theta2 = 2.0*PI*random_double();  // Random number (0, 2*PI)
	phi = 2.0*PI*random_double();     // Random number (0, 2*PI)

	if ( theta1 == theta2  && theta2 == phi ) {
	  cerr << "\n theta1, theta2 and phi are all equal.";
	} /* endif */

	aa = alpha(abs_k, theta1, phi);
	bb = beta(abs_k, theta2, phi);

	deno = abs_k * sqrt(k1*k1 + k2*k2);

	index = l + nz*(j+Ny*i);

	if (deno != 0.0) {
	  uu[index][0] = real( (aa*abs_k*k2 + bb*k1*k3)/deno );
	  uu[index][1] = imag( (aa*abs_k*k2 + bb*k1*k3)/deno );
	  
	  vv[index][0] = real( (bb*k2*k3 - aa*abs_k*k1)/deno );
	  vv[index][1] = imag( (bb*k2*k3 - aa*abs_k*k1)/deno );

	  ww[index][0] = real( -( bb*(k1*k1 + k2*k2) )/deno );
	  ww[index][1] = imag( -( bb*(k1*k1 + k2*k2) )/deno );

	} else {
	  uu[index][0] = 0.0;
	  uu[index][1] = 0.0;

	  vv[index][0] = 0.0;
	  vv[index][1] = 0.0;

	  ww[index][0] = 0.0;
	  ww[index][1] = 0.0;
	} /* endif */

 	if ( l==0  ||  l==Nz/2) {
	  // complex conjugates
	  if ( j>Ny/2  ||  ( i>Nx/2  &&  (j==0 || j==Ny/2) ) ) {
	    int index_conj = l + nz*(jconj+Ny*iconj);

	    uu[index][0] =  uu[index_conj][0];  
	    uu[index][1] = -uu[index_conj][1];

	    vv[index][0] =  vv[index_conj][0];
	    vv[index][1] = -vv[index_conj][1];

	    ww[index][0] =  ww[index_conj][0];
	    ww[index][1] = -ww[index_conj][1];

	  // real values at 8 corners
	  } else if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) {
	    uu[index][1] = 0.0;
	    vv[index][1] = 0.0;
	    ww[index][1] = 0.0; 
	  }

	} // end if

      } /* end for */
    } /* end for */
  } /* end for */
  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, uu, u, FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
 
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, vv, v, FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ww, w, FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);
  
  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int l=0; l<Nz; ++l) {
	index = l + Nz*(j+Ny*i);
	u[index] =  u[index] /*scaling_factor*/;
	v[index] =  v[index] /*scaling_factor*/;
	w[index] =  w[index] /*scaling_factor*/;
      } /* endfor */
    } /* endfor */
  } /* endfor */
  
  // Assign turbulent velocity field
  int nBlk, ix, iy, iz;
  int INl, INu, JNl, JNu, KNl, KNu;
  for (int kBlk = 0; kBlk <= Initial_Mesh.NBlk_Kdir-1; ++kBlk) {
     for (int jBlk = 0; jBlk <= Initial_Mesh.NBlk_Jdir-1; ++jBlk) {
        for (int iBlk = 0; iBlk <= Initial_Mesh.NBlk_Idir-1; ++iBlk) {
            nBlk = iBlk +
                   jBlk*Initial_Mesh.NBlk_Idir +
                   kBlk*Initial_Mesh.NBlk_Idir*Initial_Mesh.NBlk_Jdir;
            for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
               for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
                  for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
		     ix = iBlk*IPs.NCells_Turbulence_Idir+(i-Initial_Mesh.Grid_Blks[nBlk].ICl);
		     iy = jBlk*IPs.NCells_Turbulence_Jdir+(j-Initial_Mesh.Grid_Blks[nBlk].JCl);
		     iz = kBlk*IPs.NCells_Turbulence_Kdir+(k-Initial_Mesh.Grid_Blks[nBlk].KCl);
	             index = iz + 
                             iy*Nz + 
                             ix*Ny*Nz;
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].x = u[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].y = v[index];
                     Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k].z = w[index];
		     Initial_Velocity_Field.Vel_Blks[nBlk].Position[i][j][k] = Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc;	     
                     
		  } /* endfor */
	       } /* endfor */
	    } /* endfor */  
	    INl = Initial_Mesh.Grid_Blks[nBlk].INl; 
	    INu = Initial_Mesh.Grid_Blks[nBlk].INu; 
	    JNl = Initial_Mesh.Grid_Blks[nBlk].JNl; 
	    JNu = Initial_Mesh.Grid_Blks[nBlk].JNu; 
	    KNl = Initial_Mesh.Grid_Blks[nBlk].KNl; 
	    KNu = Initial_Mesh.Grid_Blks[nBlk].KNu;
	    Initial_Velocity_Field.Vel_Blks[nBlk].Node_INl_JNl_KNl = Initial_Mesh.Grid_Blks[nBlk].Node[INl][JNl][KNl].X; 
	    Initial_Velocity_Field.Vel_Blks[nBlk].Node_INu_JNu_KNu = Initial_Mesh.Grid_Blks[nBlk].Node[INu][JNu][KNu].X;

	} /* endfor */
     } /* endfor */
  } /* endfor */
  
  // Deallocations   
  fftw_free(u);
  fftw_free(v);
  fftw_free(w);
  fftw_free(uu);
  fftw_free(vv);
  fftw_free(ww);

  return (0);

}

// Write_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Write_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			       const Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                               const Grid3D_Input_Parameters &IPs) const {

    ofstream out_file;
    out_file.open("Initial_Turbulence_Fluctuations.dat", ios::out);
    if (out_file.fail()) {
      cerr<<"\n Error opening file: Initial_Turbulence_Fluctuations.dat to write" << endl;
      exit(1);
    } /* endif */
  
    out_file.setf(ios::scientific);

    for (int nBlk = 0; nBlk <= Initial_Mesh.NBlk-1; ++nBlk ) {
       for (int i = Initial_Mesh.Grid_Blks[nBlk].ICl; i <= Initial_Mesh.Grid_Blks[nBlk].ICu; ++i) {
          for (int j = Initial_Mesh.Grid_Blks[nBlk].JCl; j <= Initial_Mesh.Grid_Blks[nBlk].JCu; ++j) {
             for (int k = Initial_Mesh.Grid_Blks[nBlk].KCl; k <= Initial_Mesh.Grid_Blks[nBlk].KCu; ++k) {
	      out_file << setprecision(10)
		       << Initial_Mesh.Grid_Blks[nBlk].Cell[i][j][k].Xc << " " 
                       << Initial_Velocity_Field.Vel_Blks[nBlk].Velocity[i][j][k];
	     } /* endfor */
	  } /* endfor */
       } /* endfor */
       out_file.setf(ios::scientific);
       out_file << endl;
    } /* endfor */

    out_file.close();

}

// Read_Turbulent_Velocity_Field
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Read_Turbulent_Velocity_Field(const Grid3D_Hexa_Multi_Block_List &Initial_Mesh,
			      Turbulent_Velocity_Field_Multi_Block_List &Initial_Velocity_Field,
                              Grid3D_Input_Parameters &IPs) {

   ifstream InFile;

   // Open data file for reading
   InFile.open("Initial_Turbulence_Fluctuations.dat", ios::in); 
   // Check to see if successful
   if(InFile.fail()){ 
     cerr<<"\n Error opening file: Initial_Turbulence.dat to read" <<endl;
     exit(1); 
   } 

   bool interpolate_flag;
   double xx, yy, zz, dx, dy, dz, dd, uprime_x, uprime_y, uprime_z;

   InFile.setf(ios::skipws);
 
   InFile.unsetf(ios::skipws);
   InFile.close();

}

// Assign_Homogeneous_Turbulence_Velocity_Field
template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK *Solution_Block,
                                                  const AdaptiveBlock3D_List &LocalSolnBlockList,
                                                  const Turbulent_Velocity_Field_Multi_Block_List &Velocity_Field) {

   /* Assign initial turbulent velocity field to each solution block. */

   for (int nBlk = 0 ; nBlk <= LocalSolnBlockList.Nblk-1 ; nBlk++) {
      if (LocalSolnBlockList.Block[nBlk].used == ADAPTIVEBLOCK3D_USED) {
	 if (Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum].Allocated) {
  	    Assign_Homogeneous_Turbulence_Velocity_Field(Solution_Block[nBlk],
                                                         Velocity_Field.Vel_Blks[LocalSolnBlockList.Block[nBlk].info.gblknum]);
	 } /* endif */
      } /* endif */
   }  /* endfor */

}

template<typename HEXA_BLOCK>
void Assign_Homogeneous_Turbulence_Velocity_Field(HEXA_BLOCK &Solution_Block,
                                                  const Turbulent_Velocity_Field_Block &Velocity_Field) {

  /* Assign initial turbulent velocity field to solution block. */

  for (int i = Solution_Block.ICl ; i <= Solution_Block.ICu ; i++) {
     for (int j = Solution_Block.JCl ; j <= Solution_Block.JCu ; j++) {
        for (int k = Solution_Block.KCl ; k <= Solution_Block.KCu ; k++) {
           Solution_Block.W[i][j][k].v += Velocity_Field.Velocity[i][j][k];
 	   Solution_Block.U[i][j][k] = Solution_Block.W[i][j][k].U();
	} /* endfor */
     } /* endfor */
  } /* endfor */

}




/********************************************************
 *          Open_Turbulence_Progress_File               *
 ********************************************************/
inline int Open_Turbulence_Progress_File(ofstream &Turbulence_Progress_File,
		    		         char *File_Name,
				         const int Append_to_File) {

    int i;
    char prefix[256], extension[256], 
         turbulence_progress_file_name[256], gnuplot_file_name[256];
    char *turbulence_progress_file_name_ptr, *gnuplot_file_name_ptr;
    ofstream gnuplot_file;

    /* Determine the name of the turbulence progress file. */

    i = 0;
    while (1) {
       if (File_Name[i] == ' ' ||
           File_Name[i] == '.') break;
       prefix[i] = File_Name[i];
       i = i + 1;
       if (i > strlen(File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_turbulence_statistics");

    strcpy(extension, ".dat");
    strcpy(turbulence_progress_file_name, prefix);
    strcat(turbulence_progress_file_name, extension);

    turbulence_progress_file_name_ptr = turbulence_progress_file_name;

    /* Open the turbulence progress file. */

    if (Append_to_File) {
       Turbulence_Progress_File.open(turbulence_progress_file_name_ptr, ios::out|ios::app);
    } else {
       Turbulence_Progress_File.open(turbulence_progress_file_name_ptr, ios::out);
    } /* endif */
    if (Turbulence_Progress_File.bad()) return (1);

    /* Write the appropriate GNUPLOT command file for 
       plotting turbulence progress file information. */

    strcpy(extension, ".gplt");
    strcpy(gnuplot_file_name, prefix);
    strcat(gnuplot_file_name, extension);

    gnuplot_file_name_ptr = gnuplot_file_name;

    gnuplot_file.open(gnuplot_file_name_ptr, ios::out);
    if (gnuplot_file.fail()) return(1);

    gnuplot_file << "set title \"Turbulence parameters progress \"\n"
                 << "set xlabel \"Time \"\n"
                 << "set ylabel \"TKE/u_rms/enstrophy/Taylor/St\"\n" 
                 << "set logscale xy\n"
                 << "plot \"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"TKE\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%lf%*lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"u_rms\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%lf%*lf%*lf%*lf\" \\\n"
                 << "     title \"enstrophy\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
		 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%*lf%lf%*lf%*lf\" \\\n"
	         << "     title \"Taylor_scale\" with lines, \\\n"
                 << "\"" << turbulence_progress_file_name_ptr << "\""
                 << " using 1:2 \"%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%lf\" \\\n"
                 << "     title \"St\" with lines\n"
                 << "pause -1  \"Hit return to continue\"\n";

    gnuplot_file.close();

    /* Preparation of progress file complete.
       Return zero value. */

    return(0);

}

/********************************************************
 *          Close_Turbulence_Progress_File              *
 ********************************************************/
inline int Close_Turbulence_Progress_File(ofstream &Turbulence_Progress_File) {

    Turbulence_Progress_File.close();

    return(0);

}

/********************************************************
 *          Output_Turbulence_Progress_to_File          *
 ********************************************************/
inline int Output_Turbulence_Progress_to_File(ostream &Turbulence_Progress_File,
		 			      const int Number_of_Time_Steps,
					      const double &Time,
					      const CPUTime &CPU_Time,
					      const double &Total_Energy,
                                              const double &u_rms,
                                              const double &Total_Enstrophy,
					      const double &Taylor_scale,
					      const double &viscosity,
                                              const double &turbulent_burning_rate) {

    Turbulence_Progress_File << setprecision(6);
    Turbulence_Progress_File << Time
                             << " " << Number_of_Time_Steps
                             << " " << CPU_Time.min();
    Turbulence_Progress_File.setf(ios::scientific);
    Turbulence_Progress_File << " " << Total_Energy
                             << " " << u_rms
                             << " " << Total_Enstrophy
                             << " " << Taylor_scale
			     << " " << viscosity
			     << " " << turbulent_burning_rate
                             << "\n";                       
    Turbulence_Progress_File.unsetf(ios::scientific);
    Turbulence_Progress_File.flush();

    return(0);

}

template <typename HEXA_BLOCK>
int Longitudinal_Correlation(Octree_DataStructure &OcTree,
                             AdaptiveBlock3D_ResourceList &Global_Soln_Block_List,
			     AdaptiveBlock3D_List &LocalSolnBlockList,
			     HEXA_BLOCK *Solution_Block,
			     Grid3D_Input_Parameters &IPs,
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


#endif // _TURBULENT_VELOCITY_FIELD_INCLUDED 

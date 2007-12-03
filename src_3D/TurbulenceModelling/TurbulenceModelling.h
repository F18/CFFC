/****************** TurbulenceModelling.h **************************************
  This header file defines various classes for various standard turbulence
  models.  
********************************************************************************/

#ifndef _TURBULENCE_MODELLING_INCLUDED 
#define _TURBULENCE_MODELLING_INCLUDED

/* Include required CFFC header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

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

#ifndef _OCTREE_INCLUDED
#include "../AMR/Octree.h"
#endif // _OCTREE_INCLUDED

#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif // _FFTW_INCLUDED

// Constants
const complex<double>  I(0.0, 1.0);      // sqrt(-1.0)

class Turbulent3DWallData;
class Turbulence_Model_k_omega;

class Turbulence_Model_k_omega{
  public:
      
   static double alpha;
   static double sigma;
   static double sigma_star;
   static double beta;
   static double f_beta_const;
   static double beta_star;
   static double f_beta_star_const;
   static double EDM_const;
   static double y_sublayer;
   static double Karman_const;
   static double C_const;
   static double yplus_l;
   static double yplus_u;
   
   double omega_sublayer_BC(const double &d, const double &mu, const double &y) const;

   double f_beta(const Tensor3D &rotation_tensor, const Tensor3D &strain_rate ,
                 const double &omega) const;

   double f_betastar(const double &dkdx, const double &dkdy, const double &dkdz, 
                     const double &domegadx, const double &domegady, 
                     const double &domegadz, const double &omega) const;

   int wall_function(const double &d, const double &mu, Turbulent3DWallData &wall, 
                     double &dk, double &domega);

   int low_Reynolds_number_formulation(const double &d, const double &mu, 
                                       const int i, const int j, const int k,
                                       Turbulent3DWallData &wall,
                                       Grid3D_Hexa_Block &grid, double &domega);

   int automatic_wall_treatment(const double &d, const double &mu, 
                                const int i, const int j, const int k,
                                Turbulent3DWallData &wall, 
                                Grid3D_Hexa_Block &grid, double &dk,
                                double &domega);

   int wall_function_residual_evaluation(Turbulent3DWallData &wall, double &dUdt_rhok, double &dUdt_rhoomega);

   int low_Reynolds_number_formulation_residual_evaluaton(const int i, const int j, const int k,
                                                          Turbulent3DWallData &wall,
                                                          Grid3D_Hexa_Block &grid,
                                                          double &dUdt_rhoomega);

   int  automatic_wall_treatment_residual_evaluaton(const int i, const int j, const int k,
                                                    Turbulent3DWallData &wall,
                                                    Grid3D_Hexa_Block &grid,
                                                    double &dUdt_rhok,
                                                    double &dUdt_rhoomega);
   
};

/*!
 * Class: Turbulent3DWallData
 *
 * @brief Class defined to store and handle wall distance calculations
 * (ie, ywall, yplus, ...).
 *
 * \verbatim
 * Member functions
 *     ywall -- Distance to nearest solid wall.
 *     Xwall -- Location of the nearest solid wall.
 *     nwall -- Normal at the nearest solid wall.
 *     yplus -- Dimensionless distance to nearest wall.
 *     utau  -- Friction velocity at nearest wall.
 *     tauw  -- Shear stress at nearest wall.
 *
 * Member operators
 *      W -- a primitive solution state
 *      c -- a scalar (double)
 *
 * cout << W; (output function)
 * cin  >> W; (input function)
 * \endverbatim
 */
class Turbulent3DWallData {
  private:
  public:
   double    ywall; //!< Distance to nearest solid wall.
   Vector3D  Xwall; //!< Location of the nearest solid wall.
   Vector3D  nwall; //!< Normal at the nearest solid wall.
   double    yplus; //!< Dimensionless distance to nearest wall.
   double     utau; //!< Friction velocity at nearest wall.
   double     tauw; //!< Shear stress at nearest wall.
   int      BCwall; //!< Wall boundary condition.
   
   //! Creation constructor.
   Turbulent3DWallData(void) {
      ywall = ZERO; Xwall = Vector3D_ZERO; nwall = Vector3D_ZERO;
      yplus = ZERO; utau = ZERO; tauw = ZERO;
      BCwall = BC_WALL_VISCOUS_HEATFLUX;
   }

   //! Copy constructor.
   Turbulent3DWallData(const Turbulent3DWallData &W) {
      ywall = W.ywall; Xwall = W.Xwall; nwall = W.nwall;
      yplus = W.yplus; utau = W.utau; tauw = W.tauw;
      BCwall = W.BCwall;
   }
   
   //! Destructor.
   ~Turbulent3DWallData(void) { };
   
   //@{ @name Input-output operators.
   friend ostream &operator << (ostream &out_file, const Turbulent3DWallData &W);
   friend istream &operator >> (istream &in_file, Turbulent3DWallData &W);
   //@}

};

//! Turbulent3DWallData -- Output operator.
inline ostream &operator << (ostream &out_file, const Turbulent3DWallData &W) {
  out_file.setf(ios::scientific);
  out_file << W.ywall << W.Xwall << W.nwall << " ";
  out_file.unsetf(ios::scientific);
  out_file << W.BCwall << " ";
  out_file.setf(ios::scientific);
  out_file << W.yplus << " " << W.utau << " " << W.tauw << endl;
  out_file.unsetf(ios::scientific);
  return out_file;
}

//! Turbulent3DWallData -- Input operator.
inline istream &operator >> (istream &in_file, Turbulent3DWallData &W) {
  in_file.setf(ios::skipws);
  in_file >> W.ywall >> W.Xwall >> W.nwall >> W.BCwall >> W.yplus >> W.utau >> W.tauw;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * The next three routines determine the normal distance to the       *
 * nearest wall, the location of the nearest wall, the normal at the  *
 * nearest wall, and the boundary condition at the nearest wall for   *
 * each and every cell.                                               *
 **********************************************************************/
// template function to determine the wall distance

/**********************************************************************
 * Routine: Wall_Distance                                             *
 *                                                                    *
 * This routine determines the normal distance to the wall for each   *
 * cell of a multiblock body-fitted quadrilateral mesh.  A direct     *
 * (exhaustive) search is performed.  The wall location, the inward   *
 * normal at the wall, and the wall boundary condition type is also   *
 * found and stored.                                                  *
 *                                                                    *
 **********************************************************************/
template<class HEXA_BLOCK>
int Wall_Distance(HEXA_BLOCK &Solution_Block,
                  const Vector3D X_cell,
                  Vector3D &X_wall,
                  Vector3D &n_wall,
                  double &y_wall,
                  int &BC_wall){
   
   // Initialize y_wall.
   y_wall = 1.0e300;
   
   X_wall = Vector3D_ZERO;
   n_wall = Vector3D_ZERO;
   BC_wall = BC_WALL_VISCOUS_HEATFLUX;
   
   int i, j, k;
   
   // Check West boundary.
   for (k = Solution_Block.Grid.KCl-Solution_Block.Grid.Nghost; k <= Solution_Block.Grid.KCu+Solution_Block.Grid.Nghost; ++k) {
      for (j =  Solution_Block.Grid.JCl-Solution_Block.Grid.Nghost; j <=  Solution_Block.Grid.JCu+Solution_Block.Grid.Nghost; ++j ) {
         if (Solution_Block.Grid.BCtypeW[j][k] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeW[j][k] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeW[j][k] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeW[j][k] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceW(Solution_Block.Grid.ICl,j,k)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceW(Solution_Block.Grid.ICl,j,k)-X_cell);
               X_wall = Solution_Block.Grid.xfaceW(Solution_Block.Grid.ICl, j, k);
               n_wall = -Solution_Block.Grid.nfaceW(Solution_Block.Grid.ICl, j, k);
               BC_wall = Solution_Block.Grid.BCtypeW[j][k];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */
   
   // Check East boundary.
   for (k = Solution_Block.Grid.KCl-Solution_Block.Grid.Nghost; k <= Solution_Block.Grid.KCu+Solution_Block.Grid.Nghost; ++k) {
      for (j = Solution_Block.Grid.JCl-Solution_Block.Grid.Nghost; j <= Solution_Block.Grid.JCu+Solution_Block.Grid.Nghost; ++j) {
         if (Solution_Block.Grid.BCtypeE[j][k] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeE[j][k] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeE[j][k] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeE[j][k] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceE(Solution_Block.Grid.ICu,j,k)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceE(Solution_Block.Grid.ICu,j,k)-X_cell);
               X_wall = Solution_Block.Grid.xfaceE(Solution_Block.Grid.ICu,j,k);
               n_wall = -Solution_Block.Grid.nfaceE(Solution_Block.Grid.ICu,j,k);
               BC_wall = Solution_Block.Grid.BCtypeE[j][k];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Check South boundary.
   for (k = Solution_Block.Grid.KCl-Solution_Block.Grid.Nghost; k <= Solution_Block.Grid.KCu+Solution_Block.Grid.Nghost; ++k) {
      for (i = Solution_Block.Grid.ICl-Solution_Block.Grid.Nghost; i <= Solution_Block.Grid.ICu+Solution_Block.Grid.Nghost; ++i) {
         if (Solution_Block.Grid.BCtypeS[i][k] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeS[i][k] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeS[i][k] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeS[i][k] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceS(i, Solution_Block.Grid.JCl, k)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceS(i,Solution_Block.Grid.JCl, k)-X_cell);
               X_wall =  Solution_Block.Grid.xfaceS(i,Solution_Block.Grid.JCl, k);
               n_wall = -Solution_Block.Grid.nfaceS(i,Solution_Block.Grid.JCl, k);
               BC_wall = Solution_Block.Grid.BCtypeS[i][k];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Check North boundary.
   for (k = Solution_Block.Grid.KCl-Solution_Block.Grid.Nghost; k <=  Solution_Block.Grid.KCu+Solution_Block.Grid.Nghost; ++k) {
      for (i = Solution_Block.Grid.ICl-Solution_Block.Grid.Nghost; i <=  Solution_Block.Grid.ICu+Solution_Block.Grid.Nghost; ++i ) {
         if (Solution_Block.Grid.BCtypeN[i][k] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeN[i][k] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeN[i][k] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeN[i][k] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceN(i,Solution_Block.Grid.JCu,k)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceN(i,Solution_Block.Grid.JCu,k)-X_cell);
               X_wall = Solution_Block.Grid.xfaceN(i,Solution_Block.Grid.JCu,k);
               n_wall = -Solution_Block.Grid.nfaceN(i,Solution_Block.Grid.JCu,k);
               BC_wall = Solution_Block.Grid.BCtypeN[i][k];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Check Bottom boundary.
   for (j = Solution_Block.Grid.JCl-Solution_Block.Grid.Nghost; j <= Solution_Block.Grid.JCu+Solution_Block.Grid.Nghost; ++j) {
      for (i = Solution_Block.Grid.ICl-Solution_Block.Grid.Nghost; i <= Solution_Block.Grid.ICu+Solution_Block.Grid.Nghost; ++i) {
         if (Solution_Block.Grid.BCtypeB[i][j] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeB[i][j] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeB[i][j] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeB[i][j] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceBot(i,j,Solution_Block.Grid.KCl)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceBot(i,j,Solution_Block.Grid.KCl)-X_cell);
               X_wall = Solution_Block.Grid.xfaceBot(i,j,Solution_Block.Grid.KCl);
               n_wall = -Solution_Block.Grid.nfaceBot(i,j,Solution_Block.Grid.KCl);
               BC_wall = Solution_Block.Grid.BCtypeB[i][j];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Check Top boundary.
   for (j = Solution_Block.Grid.JCl-Solution_Block.Grid.Nghost; j <= Solution_Block.Grid.JCu+Solution_Block.Grid.Nghost; ++j) {
      for (i = Solution_Block.Grid.ICl-Solution_Block.Grid.Nghost; i <= Solution_Block.Grid.ICu+Solution_Block.Grid.Nghost; ++i) {
         if (Solution_Block.Grid.BCtypeT[i][j] == BC_WALL_VISCOUS ||
             Solution_Block.Grid.BCtypeT[i][j] == BC_NO_SLIP  ||
             Solution_Block.Grid.BCtypeT[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
             Solution_Block.Grid.BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
             Solution_Block.Grid.BCtypeT[i][j] == BC_ADIABATIC_WALL ||
             Solution_Block.Grid.BCtypeT[i][j] == BC_BURNING_SURFACE) {
            if (abs(Solution_Block.Grid.xfaceTop(i,j,Solution_Block.Grid.KCu)-X_cell) < y_wall) {
               y_wall = abs(Solution_Block.Grid.xfaceTop(i,j,Solution_Block.Grid.KCu)-X_cell);
               X_wall = Solution_Block.Grid.xfaceTop(i,j,Solution_Block.Grid.KCu);
               n_wall = -Solution_Block.Grid.nfaceTop(i,j,Solution_Block.Grid.KCu);
               BC_wall = Solution_Block.Grid.BCtypeT[i][j];
            } /* endif */
         } /* endif */
      } /* endfor */
   } /* endfor */

   // Distance to wall computed successfully.

   return 0;

}

template <class HEXA_BLOCK>
int Wall_Distance(HEXA_BLOCK *Solution_Block,
		  AdaptiveBlock3D_List &LocalSolnBlockList,
		  const Vector3D X_cell,
		  Vector3D &X_wall,
		  Vector3D &n_wall,
		  double &y_wall,
		  int &BC_wall) {
   
  int error_flag;
  double y_wall_temp;
  Vector3D X_wall_temp, n_wall_temp;
  int BC_wall_temp;

  y_wall = 1e300;
  y_wall_temp = 1e300;

  for (int nb = 0 ; nb < LocalSolnBlockList.Nblk; nb++) {
     if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK3D_USED) {
         error_flag = Wall_Distance(Solution_Block[nb],
                                    X_cell,
                                    X_wall_temp,
                                    n_wall_temp,
                                    y_wall_temp,
                                    BC_wall_temp);
       if (error_flag) return error_flag;
       if (y_wall_temp < y_wall) {
          y_wall = y_wall_temp;
          X_wall = X_wall_temp;
          n_wall = n_wall_temp;
          BC_wall = BC_wall_temp;
       } /* endif */
     } /* endif */
  } /* endfor */
  
  // Wall distance successfully found.

  return 0;

}   

/**********************************************************
 * Routine: Distance_to_Wall                              *
 *                                                        *
 * Determines the normal distance to the wall for each    *
 * cell of a multiblock body-fitted quadrilateral mesh.   *
 * A direct search is performed.                          *
 *                                                        *
 **********************************************************/
template <class HEXA_BLOCK>
int Wall_Distance(HEXA_BLOCK *Solution_Block,
                  Octree_DataStructure &Octree,
                  AdaptiveBlock3D_List &LocalSolnBlockList) {
   
  if (Solution_Block[0].Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return (0);

  int error_flag;
  int buffer_size;
  int ilow, iup, jlow, jup, klow, kup, ni, nj, nk;
  double y_wall, y_wall_temp;
  double *x_buffer;
  Vector3D x_cell, X_wall, n_wall;
  int BC_wall;

#ifdef _MPI_VERSION
  int wall_buffer_size, wall_iCPU;
  double *wall_buffer;

  wall_buffer_size = 7;
  wall_buffer = new double[wall_buffer_size];
#endif

  // Compute the distance to the nearest wall for the centroid of every cell.
  for (int iCPU = 0; iCPU < Octree.Ncpu; iCPU++) { // Loop over available processors.
    for (int iBLK = 0; iBLK < Octree.Nblk; iBLK++) { // Loop over available blocks.
      if (Octree.Blocks[iCPU][iBLK] != NULL) {
	if (Octree.Blocks[iCPU][iBLK]->block.used) {

          // Determine dimensions of block of interest.
     	  ilow = 0;
	  jlow = 0;
	  klow = 0;

	  if (LocalSolnBlockList.ThisCPU == iCPU) {
             iup = LocalSolnBlockList.Block[iBLK].info.dimen.i - 1 +
	           2*LocalSolnBlockList.Block[iBLK].info.dimen.ghost;
             jup = LocalSolnBlockList.Block[iBLK].info.dimen.j - 1 +
	           2*LocalSolnBlockList.Block[iBLK].info.dimen.ghost;
             kup = LocalSolnBlockList.Block[iBLK].info.dimen.k - 1 +
	           2*LocalSolnBlockList.Block[iBLK].info.dimen.ghost;
	  } else {
             iup = 0;
             jup = 0;
             kup = 0;
	  }
	  iup = CFFC_Maximum_MPI(iup);
	  jup = CFFC_Maximum_MPI(jup);
          kup = CFFC_Maximum_MPI(kup);

	  // Allocate buffer for storing cell centers.
	  ni = iup - ilow + 1;
	  nj = jup - jlow + 1;
          nk = kup - klow + 1;
	  x_buffer = new double[3*ni*nj*nk];
          
	  // Get coordinates of cell centers.
	  if (LocalSolnBlockList.ThisCPU == iCPU) {
             buffer_size = 0;
             for (int k = klow; k <= kup; k++) {
                for (int j = jlow; j <= jup; j++) {
                   for (int i = ilow; i <= iup; i++) {
                      x_buffer[buffer_size  ] = Solution_Block[iBLK].Grid.Cell[i][j][k].Xc.x;
                      x_buffer[buffer_size+1] = Solution_Block[iBLK].Grid.Cell[i][j][k].Xc.y;
                      x_buffer[buffer_size+2] = Solution_Block[iBLK].Grid.Cell[i][j][k].Xc.z;
                      buffer_size = buffer_size + 3;
                   } /* endfor */
                } /* endfor */
             } /* endfor */
	  } /* endif */
              
#ifdef _MPI_VERSION
	  buffer_size = 3*ni*nj*nk;
	  MPI::COMM_WORLD.Bcast(x_buffer, buffer_size, MPI::DOUBLE, iCPU);
#endif
     
	  // Calculate normal distance to the wall.
	  buffer_size = 0;
          for (int k = klow; k <= kup; k++) {
             for (int j = jlow; j <= jup; j++) {
                for (int i = ilow; i <= iup; i++) {
                   x_cell.x = x_buffer[buffer_size  ];
                   x_cell.y = x_buffer[buffer_size+1];
                   x_cell.z = x_buffer[buffer_size+2];
                   buffer_size = buffer_size + 3;
                   error_flag = Wall_Distance(Solution_Block, 
                                              LocalSolnBlockList, 
                                              x_cell, 
                                              X_wall, 
                                              n_wall, 
                                              y_wall, 
                                              BC_wall);
                   if (error_flag) return error_flag;

#ifdef _MPI_VERSION
                   // The following routine will reduce the y_wall variable
                   // to be the minimum value of y_wall on the p processors.
                   // It will also ensure that the X_wall variable is
                   // consistent with that distance.
                   y_wall_temp = CFFC_Minimum_MPI(y_wall);
                   if (fabs(y_wall - y_wall_temp) < TOLER) {
                      wall_buffer[0] = X_wall.x;
                      wall_buffer[1] = X_wall.y;
                      wall_buffer[2] = X_wall.z;
                      wall_buffer[3] = n_wall.x;
                      wall_buffer[4] = n_wall.y;
                      wall_buffer[5] = n_wall.z;
                      wall_buffer[6] = double(BC_wall);
                      wall_iCPU = LocalSolnBlockList.ThisCPU;
                   } else {
                      wall_buffer[0] = ZERO;
                      wall_buffer[1] = ZERO;
                      wall_buffer[2] = ZERO;
                      wall_buffer[3] = ZERO;
                      wall_buffer[4] = ZERO;
                      wall_buffer[5] = ZERO;
                      wall_buffer[6] = ZERO;
                      wall_iCPU = 0;
                   }
                   wall_iCPU = CFFC_Maximum_MPI(wall_iCPU);
                   MPI::COMM_WORLD.Bcast(wall_buffer, wall_buffer_size, MPI::DOUBLE, wall_iCPU);
                   X_wall = Vector3D(wall_buffer[0], wall_buffer[1], wall_buffer[2]);
                   n_wall = Vector3D(wall_buffer[3], wall_buffer[4], wall_buffer[5]);
                   y_wall = y_wall_temp;
                   BC_wall = int(wall_buffer[6]);
#endif

                   if (LocalSolnBlockList.ThisCPU == iCPU) {
                      Solution_Block[iBLK].WallData[i][j][k].ywall = y_wall;
                      Solution_Block[iBLK].WallData[i][j][k].Xwall = X_wall;
                      Solution_Block[iBLK].WallData[i][j][k].nwall = n_wall;
                      Solution_Block[iBLK].WallData[i][j][k].BCwall = BC_wall;
                   } /* endif */
                } /* endfor */
             } /* endfor */
          } /* endfor */
          
	  // Deallocate buffer for storing cell centers.
	  delete []x_buffer; x_buffer = NULL;
          
	} /* endif */
      } /* endif */
    } /* endfor */
  } /* endfor */

#ifdef _MPI_VERSION
  delete []wall_buffer; wall_buffer = NULL;
#endif

  return(0);
   
}

/*!
 * Class: RandomFieldRogallo
 *
 * \brief Class defined to generate random fluctuations using
 * Rogallo's procedure.
 *
 */
template<class SOLN_pSTATE, class SOLN_cSTATE>
class RandomFieldRogallo{
  private:
    int spectrum_flag;  //!< Turbulence kinetic energy spectrum flag.

  public:

    //! Default constructor.
    RandomFieldRogallo() : spectrum_flag(VON_KARMAN_PAO) { }

    //! Another constructor.
    RandomFieldRogallo(const int &SPECTRUM) : spectrum_flag(SPECTRUM) { }  


      double k_1(const int &n1, const double &L1) const { return /*2.0*PI*/double(n1); }///L1; }

      double k_2(const int &n2, const double &L2) const { return /*2.0*PI*/double(n2); }///L2; }

      double k_3(const int &n3, const double &L3) const { return /*2.0*PI*/double(n3); }///L3; }

    double random_double() const { return  drand48(); }


    // alpha and beta as defined by Rogallo, 1981
    complex<double> alpha(const double &abs_wave_num, 
			  const double &theta1, 
			  const double &phi) const;

    complex<double> beta(const double &abs_wave_num, 
			 const double &theta2, 
			 const double &phi) const;

    double Energy_Spectrum_Value(const double &abs_wave_num) const;

    int Generate_Velocity_Fluctuations(Grid3D_Hexa_Multi_Block &InitMeshBlks,
				       Grid3D_Input_Parameters/*<SOLN_pSTATE, SOLN_cSTATE>*/ &IPs) const;

    void Write_Initial_Turbulent_Fluctuations(Grid3D_Hexa_Multi_Block &InitMeshBlks,
					      Grid3D_Input_Parameters/*<SOLN_pSTATE, SOLN_cSTATE>*/ &IPs,
					      double *u, double *v, double *w) const;
       
};


//-----------------------------------------------------//
//         Members of RandomFieldRogallo class         //
//-----------------------------------------------------//

// alpha
template<class SOLN_pSTATE, class SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
alpha(const double &abs_wave_num, const double &theta1, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta1) * cos(phi);
}


// beta
template<class SOLN_pSTATE, class SOLN_cSTATE>
complex<double> RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
beta(const double &abs_wave_num, const double &theta2, const double &phi) const {

  double E, k;
  k = abs_wave_num;  E = Energy_Spectrum_Value(k);

  return (k == 0.0)  ?  0.0 : sqrt(E/(4.0*PI*k*k)) * exp(I*theta2) * sin(phi);
}


// Prescribed energy spectrum
template<class SOLN_pSTATE, class SOLN_cSTATE>
double RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Energy_Spectrum_Value(const double &abs_wave_num) const {

  double k, kp, kd, eps, u, Lp, EE;
  double C, A, alpha, a_s;
  int s;
 
  k = abs_wave_num;
  
  switch (spectrum_flag) {
    /*****  Lee and Reynolds  *****/
  case LEE_REYNOLDS :
    C = 20.0;  kp = 80.0;   
    EE = (k <= kp)  ?  C*k*k : C*kp*kp*pow(k/kp, -5.0/3.0);
    break;

    /*****  Laval and Nazarenko paper  *****/
  case LAVAL_NAZARENKO :
    //double EE, kp, C;
    C = 1.0;  kp = 4.0;
    EE = C*k*exp(-pow(k/kp, 2.0));
    break;

    /*****   von Karman-Pao   *****/
  case VON_KARMAN_PAO :
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp =1.0;   kd = 300.0 /*362.0*/;
    u = 0.8;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****  Haworth and Poinsot paper  *****/
  case HAWORTH_POINSOT :
    //double EE, kp, 
    u = 15.5;  Lp = TWO*PI/8.0;  // kp=4.0, u=2.5  kp=8
    kp = TWO*PI/Lp;
    EE = (32.0/3.0) * sqrt(2.0/PI)* (u*u/kp) * pow(k/kp, 4.0) * exp(-2.0*(k/kp)*(k/kp));
    break;

    /*****  Chasnov paper 1996  *****/
  case CHASNOV :  
    //double a_s, EE, kp = 4.0 /*8.0 20.0 4.0*/, u = 0.095; /*0.1 28.3 21.21  0.001*/  
    //int s = 3;
    s = 3.0;
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
  case BELL_DAY :
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

 
// Generate_Velocity_Fluctuations
template<class SOLN_pSTATE, class SOLN_cSTATE>
int RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Generate_Velocity_Fluctuations(Grid3D_Hexa_Multi_Block &InitMeshBlks,
			       Grid3D_Input_Parameters/*<SOLN_pSTATE, SOLN_cSTATE>*/ &IPs) const {



  double L1, L2, L3;
  L1 = IPs./*IP_Grid.*/Box_Length;
  L2 = IPs./*IP_Grid.*/Box_Width;
  L3 = IPs./*IP_Grid.*/Box_Height;

  int Nx, Ny, Nz;
  Nx = IPs./*IP_Grid.ICells*/NCells_Idir * IPs.NBlk_Idir;//InitMeshBlks.NBlk_Idir;//->NBI;
  Ny = IPs./*IP_Grid.JCells*/NCells_Jdir * IPs.NBlk_Jdir;//InitMeshBlks.NBlk_Jdir;//->NBJ;
  Nz = IPs./*IP_Grid.KCells*/NCells_Kdir * IPs.NBlk_Kdir;//InitMeshBlks.NBlk_Kdir;//->NBK;


  double        scaling_factor = 1.0/double(Nx*Ny*Nz);  // Scaling factor for the complex to real transform

  double        *u, *v, *w;        // Arrays to store the velocity fluctuations in physical space
  fftw_complex  *uu, *vv, *ww;     // Arrays to store the velocity fluctuations in Fourier space
  fftw_plan      physical;
    

  int nz = Nz/2+1;

  // Allocation of arrays used in the transforms
  u = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  v = (double *) malloc(Nx*Ny*Nz * sizeof(double));
  w = (double *) malloc(Nx*Ny*Nz * sizeof(double));

  uu = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  vv = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));
  ww = (fftw_complex *) fftw_malloc(Nx*Ny*nz * sizeof(fftw_complex));


  
  int seed = 3; 
  //int seed = time(NULL);   // assigns the current time to the seed
  srand48(seed);             // changes the seed for drand48()
  int iconj, jconj;          // Position of the conjugate complex for the i index
  double k1, k2, k3;         // Wave numbers
  
  double theta1, theta2, phi;
  complex<double> aa, bb;
  double deno;
  

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
	theta1 =  2.0*PI*random_double();  // Random number (0, 2*PI)
	theta2 =  2.0*PI*random_double();  // Random number (0, 2*PI)
	phi = 2.0*PI*random_double();      // Random number (0, 2*PI)

	if ( theta1 == theta2  && theta2 == phi ) {
	  cerr << "\n theta1, theta2 and phi are all equal.";
	}

	aa = alpha(abs_k, theta1, phi);
	bb = beta(abs_k, theta2, phi);

	deno = abs_k * sqrt(k1*k1 + k2*k2);

	int index = l + nz*(j+Ny*i);

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

	}

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


	// real values at 8 corners
// 	if ( (i==0 || i==Nx/2)  &&  (j==0 || j==Ny/2) ) {
	    
// 	    uu[i][j][l].im = 0.0;
// 	    vv[i][j][l].im = 0.0;
// 	    ww[i][j][l].im = 0.0;

// 	}






// 	cout << "\n (" << uu[i][j][l].re << "," << uu[i][j][l].im <<") \t"
// 	     << "(" << vv[i][j][l].re << "," << vv[i][j][l].im <<") \t"
// 	     << "(" << ww[i][j][l].re << "," << ww[i][j][l].im <<")";
	
        // cout <<"("<< k1*uu[i][j][l].re + k2*vv[i][j][l].re + k3*ww[i][j][l].re <<","
// 	     << k1*uu[i][j][l].im + k2*vv[i][j][l].im + k3*ww[i][j][l].im <<") \t";
      
        
      } // end for
      
      //cout << endl;

    } // end for
    
    //cout << endl;

  } // end for
  

  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, uu, u, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);

  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, vv, v, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);

  
  physical = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, ww, w, /*FFTW_BACKWARD,*/ FFTW_ESTIMATE);
  fftw_execute(physical); 
  fftw_destroy_plan(physical);


  
  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {
      for (int l=0; l<Nz; ++l) {
	int index = l + Nz*(j+Ny*i);
	u[index] =  u[index] /*scaling_factor*/;
	v[index] =  v[index] /*scaling_factor*/;
	w[index] =  w[index] /*scaling_factor*/;


      }
    }
  }


  Write_Initial_Turbulent_Fluctuations(InitMeshBlks, IPs, u, v, w);


  // Deallocations   
  fftw_free(u);
  fftw_free(v);
  fftw_free(w);
  fftw_free(uu);
  fftw_free(vv);
  fftw_free(ww);

  return 0;

}


// Write_Initial_Turbulent_Fluctuations
template<class SOLN_pSTATE, class SOLN_cSTATE>
void RandomFieldRogallo<SOLN_pSTATE, SOLN_cSTATE>::
Write_Initial_Turbulent_Fluctuations(Grid3D_Hexa_Multi_Block &InitMeshBlks,
				     Grid3D_Input_Parameters/*<SOLN_pSTATE, SOLN_cSTATE>*/ &IPs,
				     double *u, double *v, double *w) const {



  int ii, jj, kk;
  int Nblks_i, Nblks_j, Nblks_k;
  int ICblk, JCblk, KCblk; 
  
  Nblks_i = IPs.NBlk_Idir;//InitMeshBlks.NBlk_Idir;//->NBI;
  Nblks_j = IPs.NBlk_Jdir;//InitMeshBlks.NBlk_Jdir;//->NBJ;
  Nblks_k = IPs.NBlk_Kdir;//InitMeshBlks.NBlk_Kdir;//->NBK;

  cout << "\n\n Nblks_i: " << Nblks_i << "\t Nblks_j: " 
       << Nblks_j << "\t Nblks_k: " << Nblks_k;

  ICblk = IPs./*IP_Grid.ICells*/NCells_Idir; //  /Nblks_i;  // Cells per block in the I direction
  JCblk = IPs./*IP_Grid.JCells*/NCells_Jdir; //  /Nblks_j;  // Cells per block in the J direction
  KCblk = IPs./*IP_Grid.KCells*/NCells_Kdir; //  /Nblks_k;  // Cells per block in the J direction
  
  cout << "\n ICblk: " << ICblk << "\t JCblk: " 
       << JCblk << "\t KCblk: " << KCblk;


  int Nx, Ny, Nz;
  Nx = IPs./*IP_Grid.ICells*/NCells_Idir * Nblks_i;
  Ny = IPs./*IP_Grid.JCells*/NCells_Jdir * Nblks_j;
  Nz = IPs./*IP_Grid.KCells*/NCells_Kdir * Nblks_k;
/*   cout << "\n Nx: " << Nx << "\t Ny: " */
/*        << Ny << "\t Nz: " << Nz; */
  

  ofstream out_file;
  out_file.open("Initial_Turbulence_Fluctuations.dat", ios::out);
  if(out_file.fail()){
    cerr<<"\nError opening file: Initial_Turbulence_Fluctuations.dat to write" << endl;
    exit(1);
  }
  
    out_file.setf(ios::scientific);
    //    out_file << "VARIABLES= 'x', 'y', 'z', 'velocity-x', 'velocity-y', 'velocity-z' " << endl;

  for (int i_blk=0; i_blk<Nblks_i; ++i_blk) {
    for (int j_blk=0; j_blk<Nblks_j; ++j_blk) {
      for (int k_blk=0; k_blk<Nblks_k; ++k_blk) {
	
	int ICl = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].ICl, ICu = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].ICu;
	int JCl = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].JCl, JCu = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].JCu;
	int KCl = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].KCl, KCu = InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].KCu;

	if (ICu-ICl+1 != ICblk  ||  JCu-JCl+1 != JCblk  || KCu-KCl+1 != KCblk) {
	  cout << "\nERROR: Number of cell does not match writing Initial_Turbulence_Fluctuations.dat";
	}

	//       out_file << "ZONE   I=" << Nx/Nblks_i << "\t J=" << Ny/Nblks_j << "\t K=" << Nz/Nblks_k << endl;

	for (int i=0; i<ICblk; ++i) {
	  ii = i_blk*ICblk + i;
	  for (int j=0; j<JCblk; ++j) {
	    jj = j_blk*JCblk + j;
	    for (int k=0; k<KCblk; ++k) {
	      kk = k_blk*KCblk + k;

	      out_file << setprecision(10)
		       << InitMeshBlks.Grid_Blks[i_blk][j_blk][k_blk].Cell[i+ICl][j+JCl][k+KCl].Xc
		       << " " << u[kk + Nz*(jj+Ny*ii)]
		       << " " << v[kk + Nz*(jj+Ny*ii)]
		       << " " << w[kk + Nz*(jj+Ny*ii)] << "\n";
	    }
	  }
	}

	out_file << endl;

	
      }
    }
  }

  out_file.close();   
  

}

#endif // _TURBULENCE_MODELLING_INCLUDED 

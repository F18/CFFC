/****************** TurbulenceModelling.h **************************************
  This header file defines various classes for various standard turbulence
  models.  
********************************************************************************/

#ifndef _TURBULENCE_MODELLING_INCLUDED 
#define _TURBULENCE_MODELLING_INCLUDED

/* Include required CFFC header files. */

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif //_TENSOR3D_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "../Grid/Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_INCLUDED
#include "../AMR/AdaptiveBlock3D.h"
#endif // _ADAPTIVEBLOCK3D_INCLUDED

#ifndef _OCTREE_INCLUDED
#include "../AMR/Octree.h"
#endif // _OCTREE_INCLUDED

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

#endif // _TURBULENCE_MODELLING_INCLUDED 

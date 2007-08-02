/****************** FANS3D_ThermallyPerfectState.h **************************************
  This class defines the state variables and constructors for the euler3d thermally
  perfect mixtures
***********************************************************************/
#ifndef _TURBULENCE_MODELLING_INCLUDED 
#define _TURBULENCE_MODELLING_INCLUDED


class Turbulent3DWallData;
class Turbulence_Model_k_omega;

// Required C++ libraries
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <limits>

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif //_TENSOR3D_INCLUDED

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "../Grid/Grid3DHexaBlock.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

using namespace std;


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
                                       Grid3D_Hexa_Block *grid, double &domega);

   int automatic_wall_treatment(const double &d, const double &mu, 
                                const int i, const int j, const int k,
                                Turbulent3DWallData &wall, 
                                Grid3D_Hexa_Block *grid, double &dk,
                                double &domega);

   int wall_function_residual_evaluation(Turbulent3DWallData &wall, double &dUdt_rhok, double &dUdt_rhoomega);
   int low_Reynolds_number_formulation_residual_evaluaton(const int i, const int j, const int k,
                                                          Turbulent3DWallData &wall,
                                                          Grid3D_Hexa_Block *grid,
                                                          double &dUdt_rhoomega);

   int  automatic_wall_treatment_residual_evaluaton(const int i, const int j, const int k,
                                                    Turbulent3DWallData &wall,
                                                    Grid3D_Hexa_Block *grid,
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
   int Wall_Distance(HEXA_BLOCK *SolnBlk,
                     const Vector3D X_cell,
                     Vector3D &X_wall,
                     Vector3D &n_wall,
                     double &y_wall,
                     int &BC_wall){
   
   // Initialize y_wall.
   y_wall = 1.0e70;
   
   X_wall = Vector3D_ZERO;
   n_wall = Vector3D_ZERO;
   BC_wall = BC_WALL_VISCOUS_HEATFLUX;
   
   int i, j, k;
   
   // Check West boundary.
    
   for ( k = SolnBlk->Grid->KCl- SolnBlk->Grid->Nghost ; k <= SolnBlk->Grid->KCu+ SolnBlk->Grid->Nghost ; ++k) 
      for ( j =  SolnBlk->Grid->JCl- SolnBlk->Grid->Nghost ; j <=  SolnBlk->Grid->JCu+ SolnBlk->Grid->Nghost ; ++j ) {
         if (SolnBlk->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeW[j][k] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeW[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeW[j][k] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeW[j][k] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceW(SolnBlk->Grid->ICl,j,k)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceW(SolnBlk->Grid->ICl,j,k)-X_cell);
               X_wall = SolnBlk->Grid->xfaceW(SolnBlk->Grid->ICl, j, k);
               n_wall = -SolnBlk->Grid->nfaceW(SolnBlk->Grid->ICl, j, k);
               BC_wall = SolnBlk->Grid->BCtypeW[j][k];
            }
         }
      }
  
   
   // Check East boundary.
   for ( k =  SolnBlk->Grid->KCl- SolnBlk->Grid->Nghost ; k <= SolnBlk->Grid->KCu+ SolnBlk->Grid->Nghost ; ++k) 
      for ( j =  SolnBlk->Grid->JCl- SolnBlk->Grid->Nghost ; j <=  SolnBlk->Grid->JCu+ SolnBlk->Grid->Nghost ; ++j ) {
         
         if (SolnBlk->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeE[j][k] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeE[j][k] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeE[j][k] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeE[j][k] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceE(SolnBlk->Grid->ICu,j,k)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceE(SolnBlk->Grid->ICu,j,k)-X_cell);
               X_wall = SolnBlk->Grid->xfaceE(SolnBlk->Grid->ICu,j,k);
               n_wall = -SolnBlk->Grid->nfaceE(SolnBlk->Grid->ICu,j,k);
               BC_wall = SolnBlk->Grid->BCtypeE[j][k];
            }
         }
      }

   // Check South boundary.
   
   for ( k =  SolnBlk->Grid->KCl- SolnBlk->Grid->Nghost ; k <=  SolnBlk->Grid->KCu+ SolnBlk->Grid->Nghost ; ++k )
      for ( i =  SolnBlk->Grid->ICl- SolnBlk->Grid->Nghost ; i <=  SolnBlk->Grid->ICu+ SolnBlk->Grid->Nghost ; ++i ) {
         
         if (SolnBlk->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeS[i][k] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeS[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeS[i][k] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeS[i][k] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceS(i, SolnBlk->Grid->JCl, k)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceS(i,SolnBlk->Grid->JCl, k)-X_cell);
               X_wall =  SolnBlk->Grid->xfaceS(i,SolnBlk->Grid->JCl, k);
               n_wall = -SolnBlk->Grid->nfaceS(i,SolnBlk->Grid->JCl, k);
               BC_wall = SolnBlk->Grid->BCtypeS[i][k];
            }
         }
      }

  // Check North boundary.
   for ( k =  SolnBlk->Grid->KCl- SolnBlk->Grid->Nghost ; k <=  SolnBlk->Grid->KCu+ SolnBlk->Grid->Nghost ; ++k )
      for ( i =  SolnBlk->Grid->ICl- SolnBlk->Grid->Nghost ; i <=  SolnBlk->Grid->ICu+ SolnBlk->Grid->Nghost ; ++i ) {
         
         if (SolnBlk->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeN[i][k] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeN[i][k] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeN[i][k] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeN[i][k] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceN(i,SolnBlk->Grid->JCu,k)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceN(i,SolnBlk->Grid->JCu,k)-X_cell);
               X_wall = SolnBlk->Grid->xfaceN(i,SolnBlk->Grid->JCu,k);
               n_wall = -SolnBlk->Grid->nfaceN(i,SolnBlk->Grid->JCu,k);
               BC_wall = SolnBlk->Grid->BCtypeN[i][k];
            }
         }
      }

   // Check Bottom boundary.
   for ( j =  SolnBlk->Grid->JCl- SolnBlk->Grid->Nghost ; j <=  SolnBlk->Grid->JCu+ SolnBlk->Grid->Nghost ; ++j )
      for ( i =  SolnBlk->Grid->ICl- SolnBlk->Grid->Nghost ; i <=  SolnBlk->Grid->ICu+ SolnBlk->Grid->Nghost ; ++i ) {
         if (SolnBlk->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeB[i][j] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeB[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeB[i][j] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeB[i][j] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceBot(i,j,SolnBlk->Grid->KCl)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceBot(i,j,SolnBlk->Grid->KCl)-X_cell);
               X_wall = SolnBlk->Grid->xfaceBot(i,j,SolnBlk->Grid->KCl);
               n_wall = -SolnBlk->Grid->nfaceBot(i,j,SolnBlk->Grid->KCl);
               BC_wall = SolnBlk->Grid->BCtypeB[i][j];
            }
         }
         
      }

  // Check Top boundary.
   for ( j =  SolnBlk->Grid->JCl- SolnBlk->Grid->Nghost ; j <=  SolnBlk->Grid->JCu+ SolnBlk->Grid->Nghost ; ++j )
      for ( i =  SolnBlk->Grid->ICl- SolnBlk->Grid->Nghost ; i <=  SolnBlk->Grid->ICu+ SolnBlk->Grid->Nghost ; ++i ) {
         if (SolnBlk->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS ||
             SolnBlk->Grid->BCtypeT[i][j] == BC_NO_SLIP  ||
             SolnBlk->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS_ISOTHERMAL ||
             SolnBlk->Grid->BCtypeT[i][j] == BC_WALL_VISCOUS_HEATFLUX ||
             SolnBlk->Grid->BCtypeT[i][j] == BC_ADIABATIC_WALL ||
             SolnBlk->Grid->BCtypeT[i][j] == BC_BURNING_SURFACE) {
            if (abs(SolnBlk->Grid->xfaceTop(i,j,SolnBlk->Grid->KCu)-X_cell) < y_wall) {
               y_wall = abs(SolnBlk->Grid->xfaceTop(i,j,SolnBlk->Grid->KCu)-X_cell);
               X_wall = SolnBlk->Grid->xfaceTop(i,j,SolnBlk->Grid->KCu);
               n_wall = -SolnBlk->Grid->nfaceTop(i,j,SolnBlk->Grid->KCu);
               BC_wall = SolnBlk->Grid->BCtypeT[i][j];
            }
         }
         
      }
   

   // Distance to wall computed successfully.
   return 0;

}
 
#endif //end TURBULENCE_MODEL_INCLUDED 

/*!\file ReconstructionHelpers.h
  \brief Header file defining subroutines used in the high-order reconstruction process. */

#ifndef _RECONSTRUCTIONHELPERS_INCLUDED
#define _RECONSTRUCTIONHELPERS_INCLUDED

#include <vector>
#include "../Math/Matrix.h"
#include "../Math/Vector3D.h"
#include "../Utilities/Utilities.h"
#include "CENO_Tolerances.h"
#include "CENO_ExecutionMode.h"


/*! \brief Compute the geometric weight based on the imput distance.
 * 
 * The returned geometric weight depends on the precompilation 
 * directive _CENO_SQUARE_GEOM_WEIGHTING.
 */
inline double CENO_Geometric_Weighting(double & Distance){
  
  if (CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING){
    // calculate the weight

    if (CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING ) {
      return 1.0/(EpsilonTol::epsilon + sqr(Distance));
    }
    else {
      return 1.0/(EpsilonTol::epsilon + fabs(Distance));
    }
  } else {
    return 1.0;
  }
}


/*! \brief Compute the geometric weight based on the distance between two points.
 *
 * ControlVolumeWeight --> Storage for the geometric weight. <br>
 * DistanceBetweenCentroids -->  Distance between two points (it is POSITIVE quantity!!!). <br>
 *                               In general, this is the distance between the centroid of  <br>
 *                               the CV for which the weight is computed and the centroid of <br>
 *                               which the solution is reconstructed.                       <br>
************************************************************************************/
inline void CENO_Geometric_Weighting(double & ControlVolumeWeight, const double & DistanceBetweenCentroids){

  if (CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING){
    // calculate the weight
    ControlVolumeWeight = DistanceBetweenCentroids;

    if (CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING ) {
      ControlVolumeWeight *= ControlVolumeWeight;
    }

    ControlVolumeWeight = 1.0/(EpsilonTol::epsilon + ControlVolumeWeight);

  } else {
    ControlVolumeWeight = 1.0;
  }
}


/*! \brief Function for generating the geom coeff. for cartesian cell
******************************************************************************************/
double GeomCoeffCartesian(int p1, int p2, int p3, double deltaX, double deltaY, double deltaZ, double deltaXC, double deltaYC, double deltaZC);

/*! Set the stencil for the 3D kExact reconstruction
******************************************************************************************/
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       vector<int> & i_index, vector<int> & j_index, vector<int> & k_index));

/*! Set the stencil for the 3D kExact reconstruction
******************************************************************************************/
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       int *i_index, int *j_index, int *k_index);


// --> RR: ReconstructionHelpers.h make stencil for curved boundaries
/*! Set the stencil for the 3D kExact reconstruction used at curved boundaries
 *
 * The subroutine makes sure that enough cells are included in the stencil by adding cells
 * from the oposite directions of the boundaries.
 ******************************************************************************************/
//void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
//			       const int NorthCurvedBnd, const int SouthCurvedBnd,
//			       const int EastCurvedBnd, const int WestCurvedBnd,
//			       const int &ICl, const int &ICu, const int &JCl, const int &JCu,
//			       int & StencilDimension, 
//			       vector<int> & i_index, vector<int> & j_index);

#endif // _RECONSTRUCTIONHELPERS_INCLUDED

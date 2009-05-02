/*!\file NavierStokes2DAccuracyAssessmentMultiBlock.h
  \brief Header file defining specializations of 2D accuracy assessment multi-block class for NavierStokes solution. */

#ifndef _NAVIERSTOKES2D_ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED
#define _NAVIERSTOKES2D_ACCURACY_ASSESSMENT_2D_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h" /* Include 2D accuracy assessment
									  multi-block header file. */
#include "NavierStokes2DAccuracyAssessment.h"
#include "NavierStokes2DQuad.h"    // Include 2D Navier-Stokes quadrilateral mesh solution header file.

/*******************************************
 *          SPECIALIZATIONS                *
 ******************************************/

/*!
 * Assess the solution accuracy by calculating 
 * the lift and drag forces/coefficients.
 */
template<>
int AccuracyAssessment2D_MultiBlock::
AssessSolutionAccuracyBasedOnLiftAndDragCoefficients(NavierStokes2D_Quad_Block *SolnBlk,
						     const AdaptiveBlock2D_List &Soln_Block_List,
						     const NavierStokes2D_Input_Parameters &IP){

  try {
    
    // === Set variables for lift and drag calculation
    int NumberOfSolidBodies(NavierStokes2D_Quad_Block::GridType::BndSplineType::NumberOfSolidBodies());
    Lift.clear(); 
    Lift.assign(NumberOfSolidBodies, 0.0);
    Drag.clear(); 
    Drag.assign(NumberOfSolidBodies, 0.0);
    SolidBodyLength.clear();
    SolidBodyLength.assign(NumberOfSolidBodies, 0.0);
    FreeStreamDensity = IP.FreeStreamDensity();
    FreeStreamVelocity = IP.FreeStreamVelocity();
    double alpha;			 // angle of attack
    alpha = TWO*PI*IP.Flow_Angle/360.00; // angle of attack in radians
    double cos_alpha(cos(alpha)), sin_alpha(sin(alpha)); // cosine and sine of the angle of attach

    vector<double> Fx;
    Fx.assign(NumberOfSolidBodies, 0.0); // Net force in x-direction for each body
    vector<double> Fy;
    Fy.assign(NumberOfSolidBodies, 0.0); // Net force in y-direction for each body
    vector<double> DummySolidBodyLength; DummySolidBodyLength.assign(NumberOfSolidBodies, 0.0);

    // Counters
    int nb;

    // Test for validating points
    LiftAndDragCoeffs_Helper<NavierStokes2D_Input_Parameters> PointTest(IP);

    // === Consider special cases first ===
    if ((IP.i_Grid == GRID_FLAT_PLATE) || (IP.i_Grid == GRID_ADIABATIC_FLAT_PLATE)){
      // === Flat-plate ===
      for (nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
	if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED && SolnBlk[nb].Grid.IsThereAnySolidBoundary() ) {
	
	  if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
	    // Add skin friction forces only
	    if (SolnBlk[nb].Flow_Type != FLOWTYPE_INVISCID){
	      SolnBlk[nb].AssessAccuracy.addWallShearStressAerodynamicForcesHighOrder(Fx,Fy,
										      SolidBodyLength,
										      PointTest);
	    }
	  }
	
	} //endif
      } //endfor
      
    } else {

      // Compute the contributions to aerodynamic forces (i.e. Fx and Fy)
      // on each solid body due to each block on the current CPU.
      for (nb = 0; nb < Soln_Block_List.Nblk; ++nb) {
	if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED && SolnBlk[nb].Grid.IsThereAnySolidBoundary() ) {
	
	  if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
	    // Use the first high-order object to assess the accuracy
	    SolnBlk[nb].AssessAccuracy.addAerodynamicForcesHighOrder(Fx,Fy,SolidBodyLength);
	    // Add skin friction forces
	    if (SolnBlk[nb].Flow_Type != FLOWTYPE_INVISCID){
	      SolnBlk[nb].AssessAccuracy.addWallShearStressAerodynamicForcesHighOrder(Fx,Fy,
										      DummySolidBodyLength, //< Length already calc.
										      PointTest);
	    }
	  } else {
	    // Use the low-order reconstruction
	    SolnBlk[nb].AssessAccuracy.addAerodynamicForces(Fx,Fy,SolidBodyLength);
	    // Add skin friction forces
	    if (SolnBlk[nb].Flow_Type != FLOWTYPE_INVISCID){
	      throw runtime_error("Skin friction calculation not added for the stand numerical scheme.");
	    }
	  }

	} //endif
      } //endfor

    } // endif

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    /* =======================================================================
       =            Sum up all contributions from different blocks           =
       ======================================================================= */

    for (nb = 0; nb < NumberOfSolidBodies; ++nb){
      /* Total aerodynamic force Fx on all CPUs for each of the solid bodies */
      Fx[nb] = CFFC_Summation_MPI(Fx[nb]);
      
      /* Total aerodynamic force Fy on all CPUs for each of the solid bodies */
      Fy[nb] = CFFC_Summation_MPI(Fy[nb]);
            
      /* Total wetted surface on all CPUs for each of the solid bodies */
      SolidBodyLength[nb] = CFFC_Summation_MPI(SolidBodyLength[nb]);
    }

    // Calculate final lift and drag forces for each solid body
    for (nb = 0; nb < NumberOfSolidBodies; ++nb){
      // Lift force for the current solid body
      Lift[nb] = Fy[nb]*cos_alpha - Fx[nb]*sin_alpha;
	
      // Drag force for the current solid body
      Drag[nb] = Fy[nb]*sin_alpha + Fx[nb]*cos_alpha;
    }


    /* ==============================================
       =            Extra post-processing           =
       ============================================= */
    if (IP.i_Grid == GRID_CIRCULAR_CYLINDER){
      // Use the diameter of the cylinder to calculate the aerodynamic coefficients
      // Only one body should be present (i.e. the interior cylinder) !!!
      SolidBodyLength[0] = TWO*IP.Cylinder_Radius;
    }

    return 0;
  }
  catch (const ArgumentNullException & ){
    std::cerr << endl
	      << " ================================================================== " 
	      << endl
	      << " ERROR: The lift and drag coefficients couldn't be determined!" << endl
	      << " Please check the problem or don't require accuracy assessment!" << endl
	      << " ================================================================== " 
	      << endl;
    
    return 1;
  }

}


#endif

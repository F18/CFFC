/*!\file NavierStokes2DAccuracyAssessment.h
  \brief Header file defining specializations of 2D accuracy assessment class for NavierStokes solution. */

#ifndef _NAVIERSTOKES2D_ACCURACY_ASSESSMENT_2D_INCLUDED
#define _NAVIERSTOKES2D_ACCURACY_ASSESSMENT_2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */

/*******************************************
 *          SPECIALIZATIONS                *
 ******************************************/
template<>
inline bool LiftAndDragCoeffs_Helper<NavierStokes2D_Input_Parameters>::PointInIntegrationDomainTest(const Vector2D& Location) const{

  // Add your conditions for the particular problem here

  // === Flat-plate === 
  if (IP->i_Grid == GRID_FLAT_PLATE){
    // limit the integration up to the length of the flat-plate
    if (Location.x > IP->Plate_Length){
      return false;
    }
  }
  
  // Return 'true' by default
  return true;
}

/************************************************************************//**
 * \class _Wall_Shear_Stress_HO_Function_Wrapper_
 * \brief Adaptor for the solution block member function that calculates 
 * the wall shear stress with the high-order object.
 *
 * This adaptor is design to make the aforementioned member function of a 
 * structured solution block look like a function only of a given location and the 
 * normal at that location. \n
 * This wrapper is useful for integrating the forces due to wall shear stresses.
 ****************************************************************************/
class _Wall_Shear_Stress_HO_Function_Wrapper_{

public:
  
  /* Constructors */
  // == 2D ==
  _Wall_Shear_Stress_HO_Function_Wrapper_(NavierStokes2D_Quad_Block *SolnBlk_Ptr): SolnBlk(SolnBlk_Ptr),
										   iCell(0), jCell(0){ };

        
  //! Evaluate the high-order wall shear stress at a given location and the local normal
  double operator() (const Vector2D & GivenLocation, 
		     const Vector2D & LocalNormal);
    
  //! change cell indexes to new ones
  void NewIndexes(const int & _iCell_, const int & _jCell_){
    iCell = _iCell_;
    jCell = _jCell_;
  }

private:
  /*! Private default constructor*/
  _Wall_Shear_Stress_HO_Function_Wrapper_();

  // Local variables
  NavierStokes2D_Quad_Block *SolnBlk;	/*!< pointer to the solution block */
  int iCell, jCell;	                //!< the cell indexes in the structured solution block
};

/*!
 * Compute the wall shear stress provided by the high-order object.
 * The wall shear stress returned by this routine is oriented in 
 * the direction of the tangent, which is determined based on the LocalNormal.
 *
 * \param GivenLocation the location of interest
 * \param LocalNormal the normal direction at the given location
 */
inline double _Wall_Shear_Stress_HO_Function_Wrapper_::operator() (const Vector2D & GivenLocation, 
								   const Vector2D & LocalNormal){
  
  return SolnBlk->WallShearStress_HighOrder(iCell,jCell,
					    GivenLocation, LocalNormal);
}

/*!
 * Compute solution errors relative to the exact solution
 * set in NavierStokes2D_Quad_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \accuracy_digits the number of exact digits with which the error is sought
 * \param AccessToHighOrderVars NavierStokes2D_Quad_Block member function which returns the 
 *                              array of high-order objects in the block.
 * \param Pos  The specific high-order object which is used to compute the error.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 */
template<> inline
void AccuracyAssessment2D<NavierStokes2D_Quad_Block>::
ComputeSolutionErrorsHighOrder(const unsigned int &parameter,
			       const unsigned int &accuracy_digits,
			       const unsigned short int &Pos)
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set

    // Compute the errors
    SolnBlk->HighOrderVariable(Pos).
      ComputeSolutionErrors(wrapped_member_function_one_parameter(SolnBlk->ExactSolution(),
								  &NavierStokes2D_Quad_Block::
								  Exact_Solution_Type::SolutionForParameter,
								  parameter,
								  _dummy_param),
			    parameter,
			    accuracy_digits);
    
    // Write the final information in the designated variables
    L1() = SolnBlk->HighOrderVariable(Pos).L1();
    L2() = SolnBlk->HighOrderVariable(Pos).L2();
    LMax() = SolnBlk->HighOrderVariable(Pos).LMax();
    TotalBlockArea = SolnBlk->HighOrderVariable(Pos).BlockArea();
    CellsUsed = SolnBlk->HighOrderVariable(Pos).UsedCells();
    
  } else {
    // exact solution is not set
    throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionErrors() ERROR! There is no exact solution set!");
  }  
}

/*!
 * Compute solution errors relative to the exact solution
 * set in NavierStokes2D_Quad_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt NavierStokes2D_Quad_Block member function which returns the solution 
 *                          at a given location (PointOfInterest) using the reconstruction 
 *                          of cell (i,j) for a specified parameter.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 *
 */
template<> inline
void AccuracyAssessment2D<NavierStokes2D_Quad_Block>::
ComputeSolutionErrors(const unsigned int &parameter,
		      const unsigned int &accuracy_digits,
		      double (NavierStokes2D_Quad_Block::*ComputeSolutionAt)(const int &, const int &,
									     const Vector2D &,
									     const unsigned int &) const )
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set
    ComputeSolutionErrors(wrapped_member_function_one_parameter(SolnBlk->ExactSolution(),
								&NavierStokes2D_Quad_Block::
								Exact_Solution_Type::SolutionForParameter,
								parameter,
								_dummy_param),
			  parameter,
			  accuracy_digits,
			  ComputeSolutionAt);
  } else {
    // exact solution is not set
    throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionErrors() ERROR! There is no exact solution set!");
  }
  
}

/*!
 * Calculate the aerodynamic forces in the Cartesian x- and y-directions
 * due to skin friction using the high-order solution reconstruction.
 * The solid surfaces are detected based on the information carried by
 * the spline.
 * The forces are added to the provided variables.
 * 
 * \param Fx the aerodynamic force in x-direction
 * \param Fy the aerodynamic force in y-direction
 * \param WettedSurface the size of the surface that shows up during integration
 * \param Pos the index of the high-order variables
 */
template<>
template<> inline
void AccuracyAssessment2D<NavierStokes2D_Quad_Block>::
addWallShearStressAerodynamicForcesHighOrder(vector<double> & Fx, vector<double> & Fy, vector<double> & WettedSurface,
					     const LiftAndDragCoeffs_Helper<NavierStokes2D_Input_Parameters> & ValidateDomain,
					     const unsigned short int &Pos){

  // Visit each block boundary

  // === North Bnd
  if (SolnBlk->Grid.BndNorthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!
    SolnBlk->Grid.Integration.
      IntegratePiecewiseWallShearStressAlongBoundarySpline(NORTH,
							   _Wall_Shear_Stress_HO_Function_Wrapper_(SolnBlk),
							   Fx[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							   Fy[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							   WettedSurface[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							   ValidateDomain);
  }

  // === South Bnd
  if (SolnBlk->Grid.BndSouthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!
    SolnBlk->Grid.Integration.
      IntegratePiecewiseWallShearStressAlongBoundarySpline(SOUTH,
							   _Wall_Shear_Stress_HO_Function_Wrapper_(SolnBlk),
							   Fx[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							   Fy[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							   WettedSurface[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							   ValidateDomain);
  }

  // === East Bnd
  if (SolnBlk->Grid.BndEastSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!
    SolnBlk->Grid.Integration.
      IntegratePiecewiseWallShearStressAlongBoundarySpline(EAST,
							   _Wall_Shear_Stress_HO_Function_Wrapper_(SolnBlk),
							   Fx[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							   Fy[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							   WettedSurface[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							   ValidateDomain);
  }
  
  // === West Bnd
  if (SolnBlk->Grid.BndWestSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!
    SolnBlk->Grid.Integration.
      IntegratePiecewiseWallShearStressAlongBoundarySpline(WEST,
							   _Wall_Shear_Stress_HO_Function_Wrapper_(SolnBlk),
							   Fx[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							   Fy[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							   WettedSurface[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							   ValidateDomain);
  }
}

#endif

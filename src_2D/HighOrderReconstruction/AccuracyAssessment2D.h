/*!\file AccuracyAssessment2D.h
  \brief Header file defining 2D accuracy assessment class. */

#ifndef _ACCURACY_ASSESSMENT_2D_INCLUDED
#define _ACCURACY_ASSESSMENT_2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Vector2D.h"
#include "../Utilities/TypeDefinition.h"
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

template<typename Input_Parameters_Type>
class LiftAndDragCoeffs_Helper;

/*!
 * \class AccuracyAssessment2D
 *
 * @brief Collection of routines for assessing the solution accuracy of 2D problems
 *
 */
template<typename Quad_Soln_Block>
class AccuracyAssessment2D{
public:

  //! Constructor with solution block
  AccuracyAssessment2D(Quad_Soln_Block * AssociatedSolnBlock);

  //! Destructor
  ~AccuracyAssessment2D(void){ SolnBlk = NULL; }
 
  //! Re-associate solution block pointer
  void AssociateSolutionBlock(Quad_Soln_Block * AssociatedSolnBlock){ SolnBlk = AssociatedSolnBlock; }

  template<typename Input_Parameters_Type>
  void PrintErrorNorms(const Input_Parameters_Type & IP,
		       ostream & os);
  
  template<typename Input_Parameters_Type>
  void OutputErrorNormsTecplot(const Input_Parameters_Type & IP);
  
  void ComputeSolutionErrorsHighOrder(const unsigned int &parameter,
				      const unsigned int &accuracy_digits,
				      const unsigned short int &Pos = 0) throw(ArgumentNullException);

  void ComputeSolutionErrors(const unsigned int &parameter,
			     const unsigned int &accuracy_digits,
			     double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									  const Vector2D &,
									  const unsigned int &) const =
			     &Quad_Soln_Block::PiecewiseLinearSolutionAtLocation) throw(ArgumentNullException);
  
  template<typename Function_Object_Type>
  void ComputeSolutionErrors(Function_Object_Type FuncObj,
			     const unsigned int &parameter,
			     const unsigned int &accuracy_digits,
			     double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									  const Vector2D &,
									  const unsigned int &) const =
			     &Quad_Soln_Block::PiecewiseLinearSolutionAtLocation);

  template<typename Input_Parameters_Type>
  void ComputeSolutionEntropyErrorsHighOrder(const unsigned int &accuracy_digits,
					     const Input_Parameters_Type & IP,
					     const unsigned short int &Pos = 0) throw(ArgumentNullException);

  template<typename Input_Parameters_Type>
  void ComputeSolutionEntropyErrors(const unsigned int &accuracy_digits,
				    const Input_Parameters_Type & IP,
				    double (Quad_Soln_Block::*ComputeSolutionEntropyAt)(const int &, const int &,
											const Vector2D &,
											const unsigned int &) const =
				    &Quad_Soln_Block::SolutionEntropyAtLocation) throw(ArgumentNullException);

  /*! @brief Calculate and add to the provided variables the aerodynamic forces
    on solid bodies predicted by a high-order reconstruction */
  void addAerodynamicForcesHighOrder(vector<double> & Fx, 
				     vector<double> & Fy,
				     vector<double> & WettedSurface,
				     const unsigned short int &Pos = 0);

  /*! @brief Calculate and add to the provided variables the aerodynamic forces
    due to skin friction on solid bodies predicted by a high-order reconstruction */
  template<typename Input_Parameters_Type>
  void addWallShearStressAerodynamicForcesHighOrder(vector<double> & Fx, 
						    vector<double> & Fy,
						    vector<double> & WettedSurface,
						    const LiftAndDragCoeffs_Helper<Input_Parameters_Type> & ValidateDomain,
						    const unsigned short int &Pos = 0){
    throw runtime_error("AccuracyAssessment2D<>::addWallShearStressAerodynamicForcesHighOrder() ERROR! Specialization not implemented!!!");
  };

  /*! @brief Calculate and add to the provided variables the aerodynamic forces on 
    solid bodies predicted by a piecewise linear reconstruction */
  void addAerodynamicForces(vector<double> & Fx, 
			    vector<double> & Fy,
			    vector<double> & WettedSurface);

  //! @name Access to the error data:
  //@{
  double & L1(void) {return LNorms[0]; }     //!< return the L1 component of the error-norm vector
  const double & L1(void) const {return LNorms[0]; }     //!< return the L1 component of the error-norm vector
  double & L2(void) {return LNorms[1]; }     //!< return the L2 component of the error-norm vector
  const double & L2(void) const {return LNorms[1]; }     //!< return the L2 component of the error-norm vector
  double & LMax(void) {return LNorms[2]; }   //!< return the LMax component of the error-norm vector
  const double & LMax(void) const {return LNorms[2]; }   //!< return the LMax component of the error-norm vector
  double & BlockArea(void) {return TotalBlockArea; } //!< return the total area of the block
  const double & BlockArea(void) const {return TotalBlockArea; } //!< return the total area of the block

  double BlockL1Norm(void) { return LNorms[0]/TotalBlockArea; }	//!< return the L1 error norm for the block
  double BlockL2Norm(void) { return sqrt(LNorms[1]/TotalBlockArea); } //!< return the L2 error norm for the block
  double BlockLMaxNorm(void) { return LNorms[2]; } //!< return the LMax error norm for the block

  const unsigned int & UsedCells(void) const {return CellsUsed;}  //!< return the number of used cells for error calculation
  
  //! Return the reference entropy for any given location
  double & ReferenceEntropy(const double &x, const double &y) { return RefEntropy; }
  //@}

  //! Prepare object for a new calculation
  void ResetForNewCalculation(void){ AccuracyAssessed_Flag = false; }

private:
  Quad_Soln_Block *SolnBlk;	//!< Pointer to the solution block associated to this object

  AccuracyAssessment2D(void);	//!< Private default constructor  

  vector<double> LNorms;	//!< vector of errors/error norms
  double TotalBlockArea;	//!< the total area of the block
  bool AccuracyAssessed_Flag;   //!< internal flag used to avoid re-assessment of an already determined error
  bool Title_Error_Norms;       //!< internal flag used to ensure that the output of the Tecplot title is done only once
  bool Verbose;		        //!< internal flag controlling the screen output stream
  
  unsigned int _parameter;	//!< the state class parameter which is used for accuracy assessment
  unsigned int digits;		//!< the number of accurate digits with which the errors are evaluated

  unsigned int CellsUsed;	//!< the number of cells used for accuracy assessment

  double RefEntropy;	        //!< the reference entropy used for error calculation

  // Operating functions
  /*! @brief Compute errors for L1 norm calculation using the ComputeSolutionAt solution block member function */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL1(const int &iCell, const int &jCell,
				Function_Object_Type FuncObj,
 				double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									     const Vector2D &,
 									     const unsigned int &) const = 
 				&Quad_Soln_Block::PiecewiseLinearSolutionAtLocation);

  /*! @brief Compute errors for L2 norm calculation using the ComputeSolutionAt solution block member function */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL2(const int &iCell, const int &jCell,
				Function_Object_Type FuncObj,
 				double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
									     const Vector2D &,
 									     const unsigned int &) const = 
 				&Quad_Soln_Block::PiecewiseLinearSolutionAtLocation);
  
  template<typename Input_Parameters_Type>
  void SetVerboseFlag(const Input_Parameters_Type & IP);

};


/*!
 * Constructor the accuracy assessment object and set the link to the associated solution block
 */
template<typename Quad_Soln_Block> inline 
AccuracyAssessment2D<Quad_Soln_Block>::AccuracyAssessment2D(Quad_Soln_Block * AssociatedSolnBlock):
  SolnBlk(AssociatedSolnBlock),
  AccuracyAssessed_Flag(false), Title_Error_Norms(true), Verbose(false),
  LNorms(3,0.0),  TotalBlockArea(ZERO), _parameter(1), digits(10),
  CellsUsed(0)
{
  // 
}

/*!
 * Set the verbose flag based on the input parameters
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D<Quad_Soln_Block>::SetVerboseFlag(const Input_Parameters_Type & IP){
  if ( IP.Verbose() ){
    Verbose = true;
  } else {
    Verbose = false;
  }
}


/*!
 * Print the error norms to the "os" output stream
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D<Quad_Soln_Block>::PrintErrorNorms(const Input_Parameters_Type & IP,
							    ostream & os){
  
  // Set verbose flag
  SetVerboseFlag(IP);

  if(AccuracyAssessed_Flag == false){
    AssessSolutionAccuracy(SolnBlk,IP);
  }

  // print solution
  if (Verbose){
    os << "\n The accuracy norms based on parameter " << IP.ErrorParameter << " are:\n";
    os << " L1 = " << L1() << std::endl
       << " L2 = " << L2() << std::endl
       << " LMax = " << LMax() << std::endl;
    os.flush();
  }
}


/*!
 * Output the error norms to the output file in a tecplot format
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type>
void AccuracyAssessment2D<Quad_Soln_Block>::OutputErrorNormsTecplot(const Input_Parameters_Type & IP){
  // To be implemented
}

/*!
 * Compute solution errors relative to the exact solution
 * set in the Quad_Soln_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \accuracy_digits the number of exact digits with which the error is sought
 * \param Pos the specific high-order object which is used to compute the error.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 */
template<typename Quad_Soln_Block> inline
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrorsHighOrder(const unsigned int &parameter,
			       const unsigned int &accuracy_digits,
			       const unsigned short int &Pos)
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set

    // Compute the errors
    SolnBlk->HighOrderVariable(Pos).ComputeSolutionErrors(wrapped_member_function(SolnBlk->ExactSolution(),
										  &Quad_Soln_Block::
										  Exact_Solution_Type::Solution,
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
 * set in the Quad_Soln_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Quad_Soln_Block member function which returns the solution 
 *                          at a given location (PointOfInterest) using the reconstruction 
 *                          of cell (i,j) for a specified parameter.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 *
 */
template<typename Quad_Soln_Block> inline
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrors(const unsigned int &parameter,
		      const unsigned int &accuracy_digits,
		      double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
								   const Vector2D &,
								   const unsigned int &) const )
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set
    ComputeSolutionErrors(wrapped_member_function(SolnBlk->ExactSolution(),
						  &Quad_Soln_Block::Exact_Solution_Type::Solution,
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
 * Compute solution errors relative to the provided exact solution.
 * It is assumed that the piecewise solution representation
 * in the computational cells has been already calculated and
 * that the edges of the cells near boundaries are straight lines.
 *
 * \param FuncObj The exact solution function.
 *                It is assumed that the exact solution can take two arguments
 *                (x & y position) and returns a double.
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Quad_Soln_Block member function which returns the solution 
 *                          at a given location (PointOfInterest) using the reconstruction 
 *                          of cell (i,j) for a specified parameter.
 */
template<typename Quad_Soln_Block>
template<typename Function_Object_Type>
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrors(Function_Object_Type FuncObj,
		      const unsigned int &parameter,
		      const unsigned int &accuracy_digits,
		      double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
								   const Vector2D &,
								   const unsigned int &) const ){

  // Set the required 'parameter' and the number of accurate digits
  _parameter = parameter;
  digits = accuracy_digits;

  double CellError(0.0); // individual cell error for L1 and LMax norms
  int i,j;
  int StartI(SolnBlk->ICl),  EndI(SolnBlk->ICu),  StartJ(SolnBlk->JCl),  EndJ(SolnBlk->JCu);

  // reset errors values
  L1() = ZERO; L2() = ZERO; LMax() = ZERO;
  TotalBlockArea = ZERO;

  // Update CellsUsed
  CellsUsed = (EndI - StartI + 1)*(EndJ - StartJ + 1);

  // Assess the accuracy
  for (j = StartJ; j <= EndJ; ++j) {
    for (i = StartI; i <= EndI; ++i) {

      // Calculate the error in cell i,j for the given parameter and reconstruction
      CellError = ComputeSolutionErrorL1(i,j,
					 FuncObj,
					 ComputeSolutionAt);

      // Add current cell error contribution to L1 error norm
      L1() += CellError;

      // Add current cell error contribution to L2 error norm
      L2() += ComputeSolutionErrorL2(i,j,
				     FuncObj,
				     ComputeSolutionAt);

      // Compute block LMax Norm
      LMax() = max(LMax(),CellError/SolnBlk->Grid.Cell[i][j].A);

      // Compute the total area of the block
      TotalBlockArea += SolnBlk->Grid.Cell[i][j].A;

    }
  }

}

/*!
 * Compute solution entropy error relative to the provided reference state
 * for a high-order solution.
 *
 * \param accuracy_digits the number of accurate digits targeted by the integration routine
 * \param IP input parameters which provides the entropy reference state
 * \param Pos the specific high-order object which is used to compute the error.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type> inline
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionEntropyErrorsHighOrder(const unsigned int &accuracy_digits,
				      const Input_Parameters_Type & IP,
				      const unsigned short int &Pos)
  throw(ArgumentNullException)
{

  // Set the required number of accurate digits and the reference entropy
  digits = accuracy_digits;
  RefEntropy = IP.ReferenceEntropy();

  // Set the type of the returned value
  double _dummy_param(0.0);

  // Compute the errors
  SolnBlk->HighOrderVariable(Pos).ComputeSolutionErrors(wrapped_member_function(this,
										&AccuracyAssessment2D<Quad_Soln_Block>::
										ReferenceEntropy,
										_dummy_param),
							0, //< indicates that error is calculated based on Entropy
							accuracy_digits);
    
  // Write the final information in the designated variables
  L1() = SolnBlk->HighOrderVariable(Pos).L1();
  L2() = SolnBlk->HighOrderVariable(Pos).L2();
  LMax() = SolnBlk->HighOrderVariable(Pos).LMax();
  TotalBlockArea = SolnBlk->HighOrderVariable(Pos).BlockArea();
  CellsUsed = SolnBlk->HighOrderVariable(Pos).UsedCells();

}

/*!
 * Compute solution entropy error relative to the provided reference state.
 * It is assumed that the piecewise solution representation
 * in the computational cells has been already calculated and
 * that the edges of the cells near boundaries are straight lines.
 *
 * \param accuracy_digits the number of accurate digits targeted by the integration routine
 * \param IP input parameters which provides the entropy reference state
 * \param ComputeSolutionEntropyAt Quad_Soln_Block member function which returns the solution entropy
 *                                 at a given location (PointOfInterest) using the reconstruction of cell (i,j).
 */
template<typename Quad_Soln_Block>
template<typename Input_Parameters_Type>
void AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionEntropyErrors(const unsigned int &accuracy_digits,
			     const Input_Parameters_Type & IP,
			     double (Quad_Soln_Block::*ComputeSolutionEntropyAt)(const int &, const int &,
										 const Vector2D &,
										 const unsigned int &) const)
  throw(ArgumentNullException){

  // Set the required number of accurate digits and the reference entropy
  digits = accuracy_digits;
  RefEntropy = IP.ReferenceEntropy();

  double CellError(0.0); // individual cell error for L1 and LMax norms
  int i,j;
  int StartI(SolnBlk->ICl),  EndI(SolnBlk->ICu),  StartJ(SolnBlk->JCl),  EndJ(SolnBlk->JCu);

  // reset errors values
  L1() = ZERO; L2() = ZERO; LMax() = ZERO;
  TotalBlockArea = ZERO;

  // Update CellsUsed
  CellsUsed = (EndI - StartI + 1)*(EndJ - StartJ + 1);


  // Assess the accuracy
  for (j = StartJ; j <= EndJ; ++j) {
    for (i = StartI; i <= EndI; ++i) {

      // Calculate the error in cell i,j for the given parameter and reconstruction
      CellError = ComputeSolutionErrorL1(i,j,
					 wrapped_member_function(this,
								 &AccuracyAssessment2D<Quad_Soln_Block>::ReferenceEntropy,
								 CellError),
					 ComputeSolutionEntropyAt);

      // Add current cell error contribution to L1 error norm
      L1() += CellError;

      // Add current cell error contribution to L2 error norm
      L2() += ComputeSolutionErrorL2(i,j,
				     wrapped_member_function(this,
							     &AccuracyAssessment2D<Quad_Soln_Block>::ReferenceEntropy,
							     CellError),
				     ComputeSolutionEntropyAt);

      // Compute block LMax Norm
      LMax() = max(LMax(),CellError/SolnBlk->Grid.Cell[i][j].A);

      // Compute the total area of the block
      TotalBlockArea += SolnBlk->Grid.Cell[i][j].A;

    }
  }
}


/*!
 * Compute the solution error used for calculation of L1 and LMax norms
 * based on the provided exact solution and the piecewise solution 
 * representation in the computational cell (iCell,jCell).
 * The edges of the quadrilateral cell (iCell,jCell) are treated as straight lines.
 * 
 * \param FuncObj The exact solution function.
 *                It is assumed that the exact solution can take two arguments
 *                (x & y position) and returns a double.
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Quad_Soln_Block member function which returns the numerical 
 *                          solution at a given location (PointOfInterest) for cell (iCell,jCell).
 */
template<typename Quad_Soln_Block>
template<typename Function_Object_Type>
double AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrorL1(const int &iCell, const int &jCell,
		       Function_Object_Type FuncObj,
		       double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
								    const Vector2D &,
								    const unsigned int &) const ){
  // Set the type of the returned value
  double _dummy_param(0.0);
  Vector2D _dummy_Position(0.0);

  // Call the integration function
  return ( SolnBlk->Grid.
	   Integration.IntegrateFunctionOverCell(iCell,jCell,
						 error_function(FuncObj,
								wrapped_member_function_one_parameter(SolnBlk,
												      ComputeSolutionAt,
												      _dummy_Position,
												      iCell, jCell,
												      _parameter,
												      _dummy_param),
								_dummy_param),
						 digits,_dummy_param) );
  
}


/*!
 * Compute the solution error used for calculation of L2 norm
 * based on the provided exact solution and the piecewise solution 
 * representation in the computational cell (iCell,jCell).
 * The edges of the quadrilateral cell (iCell,jCell) are treated as straight lines.
 * 
 * \param FuncObj The exact solution function.
 *                It is assumed that the exact solution can take two arguments
 *                (x & y position) and returns a double.
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Quad_Soln_Block member function which returns the numerical 
 *                          solution at a given location (PointOfInterest) for cell (iCell,jCell).
 */
template<typename Quad_Soln_Block>
template<typename Function_Object_Type>
double AccuracyAssessment2D<Quad_Soln_Block>::
ComputeSolutionErrorL2(const int &iCell, const int &jCell,
		       Function_Object_Type FuncObj,
		       double (Quad_Soln_Block::*ComputeSolutionAt)(const int &, const int &,
								    const Vector2D &,
								    const unsigned int &) const ){

  // Set the type of the returned value
  double _dummy_param(0.0);
  Vector2D _dummy_Position(0.0);

  // Call the integration function
  return ( SolnBlk->Grid.
	   Integration.IntegrateFunctionOverCell(iCell,jCell,
						 square_error_function(FuncObj,
								       wrapped_member_function_one_parameter(SolnBlk,
													     ComputeSolutionAt,
													     _dummy_Position,
													     iCell, jCell,
													     _parameter,
													     _dummy_param),
								       _dummy_param),
						 digits,_dummy_param) );
}

/*!
 * Calculate the aerodynamic forces in the Cartesian x- and y-directions
 * using the high-order solution reconstruction.
 * The solid surfaces are detected based on the information carried by
 * the spline.
 * The forces are added to the provided variables.
 * Aerodynamic forces are obtained by integrating the product of between 
 * the pressure distrubtion and the local normal direction along the 
 * contour of interest.
 * 
 * \param Fx the aerodynamic force in x-direction
 * \param Fy the aerodynamic force in y-direction
 * \param WettedSurface the size of the surface that shows up during integration
 * \param Pos the index of the high-order variables
 */
template<typename Quad_Soln_Block>
void AccuracyAssessment2D<Quad_Soln_Block>::
addAerodynamicForcesHighOrder(vector<double> & Fx, vector<double> & Fy, vector<double> & WettedSurface,
			      const unsigned short int &Pos){
  
  double _dummy_param(0);

  // Define high-order data type
  typedef typename Quad_Soln_Block::HighOrderType HighOrderType;

  // Visit each block boundary

  // === North Bnd
  if (SolnBlk->Grid.BndNorthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(NORTH,
							      wrapped_soln_block_member_function(&SolnBlk->HighOrderVariable(Pos),
												 &HighOrderType::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1]);
  }

  // === South Bnd
  if (SolnBlk->Grid.BndSouthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(SOUTH,
							      wrapped_soln_block_member_function(&SolnBlk->HighOrderVariable(Pos),
												 &HighOrderType::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1]);
  }

  // === East Bnd
  if (SolnBlk->Grid.BndEastSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(EAST,
							      wrapped_soln_block_member_function(&SolnBlk->HighOrderVariable(Pos),
												 &HighOrderType::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndEastSpline.getBodyID() - 1]);
  }
  
  // === West Bnd
  if (SolnBlk->Grid.BndWestSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(WEST,
							      wrapped_soln_block_member_function(&SolnBlk->HighOrderVariable(Pos),
												 &HighOrderType::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndWestSpline.getBodyID() - 1]);
  }
}


/*!
 * Calculate the aerodynamic forces in the Cartesian x- and y-directions
 * using the piecewise linear solution reconstruction.
 * 
 * \param Fx the aerodynamic force in x-direction
 * \param Fy the aerodynamic force in y-direction
 * \param WettedSurface the size of the surface that shows up during integration
 */
template<typename Quad_Soln_Block>
void AccuracyAssessment2D<Quad_Soln_Block>::addAerodynamicForces(vector<double> & Fx, 
								 vector<double> & Fy,
								 vector<double> & WettedSurface){

  double _dummy_param(0);

  // Visit each block boundary

  // === North Bnd
  if (SolnBlk->Grid.BndNorthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(NORTH,
							      wrapped_soln_block_member_function(SolnBlk,
												 &Quad_Soln_Block::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndNorthSpline.getBodyID() - 1]);
  }

  // === South Bnd
  if (SolnBlk->Grid.BndSouthSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(SOUTH,
							      wrapped_soln_block_member_function(SolnBlk,
												 &Quad_Soln_Block::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndSouthSpline.getBodyID() - 1]);
  }

  // === East Bnd
  if (SolnBlk->Grid.BndEastSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(EAST,
							      wrapped_soln_block_member_function(SolnBlk,
												 &Quad_Soln_Block::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndEastSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndEastSpline.getBodyID() - 1]);
  }
  
  // === West Bnd
  if (SolnBlk->Grid.BndWestSpline.IsSolidBoundary()){
    // Pass dummy cell indexes (0,0) to the wrapper.
    // They will be changed correctly by the integration routine!!!

    SolnBlk->Grid.Integration.
      IntegratePiecewiseFunctionProjectionAlongBoundarySpline(WEST,
							      wrapped_soln_block_member_function(SolnBlk,
												 &Quad_Soln_Block::
												 SolutionPressureAtCoordinates,
												 0, 0, 
												 _dummy_param),
							      Fx[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							      Fy[SolnBlk->Grid.BndWestSpline.getBodyID() - 1],
							      WettedSurface[SolnBlk->Grid.BndWestSpline.getBodyID() - 1]);
  }
}


/*!
 * Class which provides useful routines to  
 * the calculation of lift and drag coefficients.
 */
template<typename Input_Parameters_Type>
class LiftAndDragCoeffs_Helper{

public:
  LiftAndDragCoeffs_Helper(const Input_Parameters_Type & IP_Ptr): IP(&IP_Ptr){ };

  //! @brief This routine tests whether the given location belongs to the integration domain.
  bool PointInIntegrationDomainTest(const Vector2D& Location) const {
    throw runtime_error("LiftAndDragCoeffs_Helper::PointInIntegrationDomainTest() ERROR! Specialize this routine!");
  }

private:
  //! Private default constructor
  LiftAndDragCoeffs_Helper(void){};

  const Input_Parameters_Type * IP; //! < pointer to input parameters
};

#endif

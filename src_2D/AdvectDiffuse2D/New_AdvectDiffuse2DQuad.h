/*!\file New_AdvectDiffuse2DQuad.h
  \brief Temporary header file defining 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

#ifndef _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED
#define _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "New_AdvectDiffuse2DState.h"  /* Include 2D advection diffusion equation solution state header file */
#include "../Grid/Cell2D.h"            /* Include 2D cell header file */
#include "../Grid/Grid2DQuad.h"        /* Include 2D quadrilateral multiblock grid header file */
#include "../AMR/QuadTree.h"           /* Include quadtree header file */
#include "../AMR/AMR.h"                /* Include AMR header file */
#include "AdvectDiffuse2DInput.h"      /* Include 2D advection diffusion equation header file */
#include "../ICEM/ICEMCFD.h"           /* Include ICEMCFD header file. */
#include "../System/System_Linux.h"    /* Include System Linux header file. */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D    3

/*!
 * Class: AdvectDiffuse2D_Quad_Block_New
 *
 * @brief Class definition of the 2D advection-diffusion solution blocks.
 *
 * \verbatim
 * Member functions
 *       U      -- Return solution variable for the block (cell average).
 *    Grid      -- Return the solution block quadrilateral grid.
 *      dt      -- Return local time step for the solution block.
 *    dudt      -- Return the solution residuals.
 *    dUdx      -- Return the unlimited solution gradients (x-direction)
 *                 for block.
 *    dUdy      -- Return the unlimited solution gradients (y-direction)
 *                 for block.
 *     phi      -- Return the solution slope limiter.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in
 *                 the i-direction (zeta-direction).
 *     ICl      -- Return lower index for cells in 
 *                 the i-direction (zeta-direction).
 *     ICu      -- Return upper index for cells in 
 *                 the i-direction (zeta-direction).
 *     NCj      -- Return number of cells in 
 *                 the j-direction (eta-direction).
 *     JCl      -- Return lower index for cells in 
 *                 the j-direction (eta-direction).
 *     JCu      -- Return upper index for cells in 
 *                 the j-direction (eta-direction).
 *  Nghost      -- Return number of ghost (halo or 
 *                 overlap) cells.
 * Axisymmetric -- Return axisymmetric geometry
 *                 indicator (=1 for axisymmetric,
 *                            =0 for planar geometry).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                 (=1 for limiter freezing on,
 *                 (=0 for limiter freezing off).
 *     UoN      -- Return array of reference states for the application
 *                 of boundary conditions at the north boundary of the
 *                 solution block.
 *     UoS      -- Return array of reference states for the application
 *                 of boundary conditions at the south boundary of the
 *                 solution block.
 *     UoE      -- Return array of reference states for the application
 *                 of boundary conditions at the east boundary of the
 *                 solution block.
 *     UoW      -- Return array of reference states for the application
 *                 of boundary conditions at the west boundary of the
 *                 solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Un      -- Return solution state at the specified node.
 *    UnNW      -- Return solution state at the north-west node.
 *    UnNE      -- Return solution state at the north-east node.
 *    UnSW      -- Return solution state at the south-west node.
 *    UnSE      -- Return solution state at the south-east node.
 *      un      -- Return solution value at the specified node.
 *    unNW      -- Return solution value at the north-west node.
 *    unNE      -- Return solution value at the north-east node.
 *    unSW      -- Return solution value at the south-west node.
 *    unSE      -- Return solution value at the south-east node.
 * evaluate_limiters -- Set flag to evaluate limiters.
 * freeze_limiters -- Set flag to freeze limiters.
 * Member functions required for message passing.
 *    NumVar    -- Returns number of solution variables in primitive
 *                 and conserved solution state vectors.
 * LoadSendBuffer -- Loads send buffer.
 * LoadSendBuffer_F2C -- Loads send buffer for fine to coarse block
 *                       messages.
 * UnloadReceiveBuffer -- Unloads receive buffer.
 * UnloadReceiveBuffer_F2C -- Unloads receive buffer for fine to
 *                            coarse block messages.
 * SubcellReconstruction -- Performs subcell solution reconstruction
 *                          used in adaptive mesh refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending conservative
 *                            flux corrections from fine to coarse
 *                            solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending 
 *                              conservative flux corrections from fine 
 *                              to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D advection diffusion equation solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_Quad_Block_New{
public:

  //! @name Solution state arrays:
  //@{
  AdvectDiffuse2D_State_New    **U; //!< Solution state.
  //@}

  //! @name Grid block information:
  //@{
  int                      NCi, //!< Total number of i-direction cells.
                           ICl, //!< First i-direction non-ghost cell counter.
                           ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
                           JCl, //!< First j-direction non-ghost cell counter.
                           JCu; //!< Final j-direction non-ghost cell counter.
  int                   Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block       Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //! @name Residual and time-stepping arrays:
  //@{
  double                         **dt; //!< Local time step.
  AdvectDiffuse2D_State_New   ***dUdt; //!< Solution residual.
  AdvectDiffuse2D_State_New      **Uo; //!< Initial solution.
  static int residual_variable;        //!< Static integer that indicates which variable is used for residual calculations.
  static int Number_of_Residual_Norms; //!< How many Residual norms to plot?
  //@}

  //! @name Solution gradient arrays:
  //@{
  AdvectDiffuse2D_State_New    **dUdx; //!< Unlimited solution gradient (x-direction).
  AdvectDiffuse2D_State_New    **dUdy; //!< Unlimited solution gradient (y-direction).
  AdvectDiffuse2D_State_New     **phi; //!< Solution slope limiter.
  //@}

  //! @name Boundary solution flux arrays:
  //@{
  AdvectDiffuse2D_State_New   *FluxN, //!< North boundary solution flux.
                              *FluxS, //!< South boundary solution flux.
                              *FluxE, //!< East boundary solution flux.
                              *FluxW; //!< West boundary solution flux.
  //@}

  //! @name Problem indicator flags:
  //@{
  static int      Axisymmetric; //!< Axisymmetric geometry indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static int         Flow_Type; //!< Flow type flag ( required by AMR ).
  //@}

  //! @name Boundary condtion reference states:
  //@{
  AdvectDiffuse2D_State_New   *UoN, //!< Boundary condition reference states for north boundary.
                              *UoS, //!< Boundary condition reference states for south boundary.
                              *UoE, //!< Boundary condition reference states for east boundary.
                              *UoW; //!< Boundary condition reference states for west boundary.
  //@}

  //! @name Accuracy assessment data:
  //@{
  static AdvectDiffuse2D_ExactSolutions *ExactSoln;          //!< Pointer to the exact solution

  /*! Variable that provides access to accuracy assessment subroutines */
  AccuracyAssessment2D<AdvectDiffuse2D_Quad_Block_New> AssessAccuracy;
  //@}
	      
  //! @name Creation, copy, and assignment constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_Quad_Block_New(void);

  /* Destructor. */
  ~AdvectDiffuse2D_Quad_Block_New(void){ deallocate(); }
  //@}

  //! @name Allocate and deallocate functions.
  //@{
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //! @name Bilinear interplation (Zingg & Yarrow).
  //@{
  //! Return solution state at specified node.
  AdvectDiffuse2D_State_New Un(const int ii, const int jj);

  AdvectDiffuse2D_State_New UnNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  AdvectDiffuse2D_State_New UnNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  AdvectDiffuse2D_State_New UnSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  AdvectDiffuse2D_State_New UnSW(const int ii, const int jj); //!< Return solution state at cell SW node.

  //! Return solution state at specified node.
  double un(const int ii, const int jj);

  double unNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  double unNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  double unSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  double unSW(const int ii, const int jj); //!< Return solution state at cell SW node.
  //@}

  //! @name Evaluate diffusive flux for the cell.
  //@{
  //@}

  //! @name Source term member functions
  //@{
  AdvectDiffuse2D_State_New SourceTerm(const int & ii, const int & jj) const;
  AdvectDiffuse2D_State_New AxisymmetricSourceTerm(const int & ii, const int & jj) const;
  //! Determine the stability limit imposed by the source term
  double SourceTermStabilityLimit(const int & ii, const int & jj) const;
  //@}

  //! @name Evaluate velocity and diffusion coefficient at different locations
  //@{
  Vector2D VelocityAtLocation(const Vector2D &PointOfInterest) const;        //!< Evaluate velocity at the PointOfInterest location
  Vector2D VelocityAtCellCentroid(const int & ii, const int & jj) const;     //!< Calculate velocity at the cell centroid (ii,jj)
  double DiffusionCoeffAtCellCentroid(const int & ii, const int & jj) const; //!< Calculate diffusion coeff. at the cell centroid
  //@}

  //! @name Member functions for limiter freezing.
  //@{
  void evaluate_limiters(void);
  void freeze_limiters(void);
  //@}

  //! @name Member functions to compute the piecewise linear solution at a particular location
  //@{
  AdvectDiffuse2D_State_New PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
							    const double &DeltaXToCentroid,
							    const double &DeltaYToCentroid) const;
  double PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
					 const double &DeltaXToCentroid,
					 const double &DeltaYToCentroid,
					 const unsigned int &parameter) const;
  AdvectDiffuse2D_State_New PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
							      const Vector2D &CalculationPoint) const;
  double PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
					   const Vector2D &CalculationPoint,
					   const unsigned int &parameter) const;
  //@}

  //! @name Member functions to compute the inviscid flux state at a boundary interface
  //@{
  void InviscidFluxStateAtBoundaryInterface(const int &BOUNDARY,
					    const int &ii, const int &jj,
					    AdvectDiffuse2D_State_New &Ul,
					    AdvectDiffuse2D_State_New &Ur);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const AdvectDiffuse2D_Quad_Block_New &Soln);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Quad_Block_New &Soln);
  //@}

  //! @name Member functions required for message passing.
  //@{
  //! Number of solution state variables.
  int NumVar(void);
  //! Load send message passing buffer.
  int LoadSendBuffer(double *buffer,
		     int &buffer_count,
		     const int buffer_size,
		     const int i_min, 
		     const int i_max,
		     const int i_inc,
		     const int j_min, 
		     const int j_max,
		     const int j_inc);
  //! Load F2C send message passing buffer.
  int LoadSendBuffer_F2C(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min, 
			 const int i_max,
			 const int i_inc,
			 const int j_min, 
			 const int j_max,
			 const int j_inc);
  //! Load C2F send message passing buffer.
  int LoadSendBuffer_C2F(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min, 
			 const int i_max,
			 const int i_inc,
			 const int j_min, 
			 const int j_max,
			 const int j_inc,
			 const int face,
			 const int sector);
  //! Unload receive message passing buffer.
  int UnloadReceiveBuffer(double *buffer,
			  int &buffer_count,
			  const int buffer_size,
			  const int i_min, 
			  const int i_max,
			  const int i_inc,
			  const int j_min, 
			  const int j_max,
			  const int j_inc);
  //! Unload F2C receive message passing buffer.
  int UnloadReceiveBuffer_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Unload C2F receive message passing buffer.
  int UnloadReceiveBuffer_C2F(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Subcell solution reconstruction within given computational cell.
  void SubcellReconstruction(const int i,
			     const int j,
			     const int Limiter);
  //! Load F2C conservative flux message passing buffer.
  int LoadSendBuffer_Flux_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min, 
			      const int i_max,
			      const int i_inc,
			      const int j_min, 
			      const int j_max,
			      const int j_inc);
  //! Unload F2C conservative flux message passing buffer.
  int UnloadReceiveBuffer_Flux_F2C(double *buffer,
				   int &buffer_count,
				   const int buffer_size,
				   const int i_min, 
				   const int i_max,
				   const int i_inc,
				   const int j_min, 
				   const int j_max,
				   const int j_inc);
  //@}

private:
  AdvectDiffuse2D_Quad_Block_New(const AdvectDiffuse2D_Quad_Block_New &Soln); //!< Private copy constructor
  AdvectDiffuse2D_Quad_Block_New operator = (const AdvectDiffuse2D_Quad_Block_New &Soln);   //!< Private assignment operator
};


/*!
 * Class: SourceTermFunctionalWithPiecewiseLinear
 *
 * @brief Definition of the source term functional
 *        with piecewise linear reconstruction
 *
 * An object of this class calculates the non-linear source 
 * term at a location specified by the Cartesian coordinates (x,y)
 * using information from the cell (ii,jj).
 * The piecewise linear reconstruction of the cell (ii,jj) is used 
 * for providing the value of the solution at the given location.
 *
 */
class SourceTermFunctionalWithPiecewiseLinear{
public:
  /*! Constructor with the cell indexes
   *  \param _ii_      the i-index of the cell that provides the information
   *  \param _jj_      the j-index of the cell that provides the information
   *  \param _SolnBlk_ the pointer to the solution block
   */
  SourceTermFunctionalWithPiecewiseLinear(const int &_ii_,
					  const int& _jj_,
					  const AdvectDiffuse2D_Quad_Block_New *_SolnBlk_): ii(_ii_),
											    jj(_jj_),
											    SolnBlk(_SolnBlk_){ };
  
  //! Return the pointwise value of the source term at a specified location (x,y)
  double operator()(const double &x, const double &y){
    return SolnBlk->U[ii][jj].s(x,y,SolnBlk->PiecewiseLinearSolutionAtLocation(ii,jj,Vector2D(x,y),1));
  }
  
private:
  int ii,jj;	  //!< cell indexes
  const AdvectDiffuse2D_Quad_Block_New * SolnBlk; //!< pointer to the solution block
  
  //! private default cstr
  SourceTermFunctionalWithPiecewiseLinear();
  //! private assignment operator
  SourceTermFunctionalWithPiecewiseLinear operator=(const SourceTermFunctionalWithPiecewiseLinear &);
};


/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::Un?? -- Get cell node solution states.     *
 **************************************************************************/
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::un?? -- Get cell node solution values.     *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::unNW(const int ii, const int jj) {
  return (un(ii, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unNE(const int ii, const int jj) {
  return (un(ii+1, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSE(const int ii, const int jj) {
  return (un(ii+1, jj));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSW(const int ii, const int jj) {
  return (un(ii, jj));
}

/*********************//**
 * Node solution value. 
 ***********************/
inline double AdvectDiffuse2D_Quad_Block_New::un(const int ii, const int jj) {
  return Un(ii,jj)[1];
}

/***************************************************//**
 * Calculate velocity at the PointOfInterest location
 ****************************************************/
inline Vector2D AdvectDiffuse2D_Quad_Block_New::VelocityAtLocation(const Vector2D &PointOfInterest) const {
  return (*U)->V(PointOfInterest.x,PointOfInterest.y);
}

/***************************************************//**
 * Calculate velocity at the centroid of cell (ii,jj).
 ****************************************************/
inline Vector2D AdvectDiffuse2D_Quad_Block_New::VelocityAtCellCentroid(const int & ii, const int & jj) const {
  return U[ii][jj].V(Grid.XCellCentroid(ii,jj), Grid.YCellCentroid(ii,jj));
}

/**************************************************************//**
 * Calculate diffusion coefficient at the centroid of cell (ii,jj).
 * Use the solution stored in the state class as
 * parameter for the diffusion coefficient field.
 ******************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::DiffusionCoeffAtCellCentroid(const int & ii, const int & jj) const {
  return U[ii][jj].k(Grid.XCellCentroid(ii,jj), Grid.YCellCentroid(ii,jj), U[ii][jj].u);
}

/**************************************************************//**
 * Determine the stability limit imposed by the source term
 * Use the solution stored in the state class as the parameter 
 * at the centroid of cell (ii,jj) for the non-linear source term field
 ******************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::SourceTermStabilityLimit(const int & ii, const int & jj) const {
  return U[ii][jj].SourceTerm->getStabilityLimit(Grid.XCellCentroid(ii,jj),
						 Grid.YCellCentroid(ii,jj),
						 U[ii][jj].u);
}

/**************************************************************//**
 * Determine the axisymmetric source term corresponding to cell (ii,jj).
 * This routine returns actually the integral of the axisymmetric source
 * term over the domain of cell (ii,jj) divided by the area.
 ******************************************************************/
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::AxisymmetricSourceTerm(const int & ii,
											const int & jj) const{
  return AdvectDiffuse2D_State_New(U[ii][jj].s_axi(Grid.Cell[ii][jj].Xc));
}


/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::evaluate_limiters -- Set flag to evaluate limiters.
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::freeze_limiters -- Set flag to freeze limiters.    *
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::NumVar -- Returns number of state variables.    *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::NumVar(void) {
  return (int(NUM_VAR_ADVECTDIFFUSE2D));
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 */
inline AdvectDiffuse2D_State_New 
AdvectDiffuse2D_Quad_Block_New::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
								const double &DeltaXToCentroid,
								const double &DeltaYToCentroid) const{
  return U[ii][jj] + (phi[ii][jj]^(dUdx[ii][jj]*DeltaXToCentroid + dUdy[ii][jj]*DeltaYToCentroid) );
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * for a particular state parameter at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 * \param parameter the parameter for which the solution is computed
 */
inline double AdvectDiffuse2D_Quad_Block_New::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
									      const double &DeltaXToCentroid,
									      const double &DeltaYToCentroid,
									      const unsigned int &parameter) const{
  return PiecewiseLinearSolutionForDelta(ii,jj,DeltaXToCentroid,DeltaYToCentroid)[parameter];
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline AdvectDiffuse2D_State_New 
AdvectDiffuse2D_Quad_Block_New::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
								  const Vector2D &CalculationPoint) const{
  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	
  return PiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y);
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * for a particular state parameter at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 * \param parameter the parameter for which the solution is computed
 */
inline double AdvectDiffuse2D_Quad_Block_New::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
										const Vector2D &CalculationPoint,
										const unsigned int &parameter) const{
  return PiecewiseLinearSolutionAtLocation(ii,jj,CalculationPoint)[parameter];
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Single Block External Subroutines.       *
 **************************************************************************/
extern void Write_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk1,
		                const AdvectDiffuse2D_Quad_Block_New &SolnBlk2);

extern int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Fine,
				  AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Coarse,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SE,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NE);

extern void Output_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
			   AdvectDiffuse2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				 AdvectDiffuse2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void ICs(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
 	        const AdvectDiffuse2D_Input_Parameters &IP,
                AdvectDiffuse2D_State_New *Wo);

extern void Set_Boundary_Reference_State(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					 const AdvectDiffuse2D_Input_Parameters &IP);

extern void BCs(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
		const AdvectDiffuse2D_Input_Parameters &IP);

extern double CFL(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  AdvectDiffuse2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                            const int i_stage,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                                       const int i_stage,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters);


/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines.     *
 **************************************************************************/

extern AdvectDiffuse2D_Quad_Block_New* Allocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
						AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern AdvectDiffuse2D_Quad_Block_New* Deallocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
						  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void ICs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
		AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern void BCs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern double CFL(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr, 
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l1_norm);

extern void L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l2_norm);

extern void Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      double *max_norm);

extern void Evaluate_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

/********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines for Mesh.  *
 ********************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Solvers.                                 *
 **************************************************************************/

extern int New_AdvectDiffuse2DQuadSolver(char *Input_File_Name_ptr,
					 int batch_flag);

#endif /* _ADVECTDIFFUSE2D_QUAD_INCLUDED  */

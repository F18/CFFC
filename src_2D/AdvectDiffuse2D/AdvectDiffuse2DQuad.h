/*!\file AdvectDiffuse2DQuad.h
  \brief Header file defining 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#define _ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DState.h"      /* Include 2D advection diffusion equation solution state header file */
#include "../Grid/Cell2D.h"            /* Include 2D cell header file */
#include "../Grid/Grid2DQuad.h"        /* Include 2D quadrilateral multiblock grid header file */
#include "../Grid/HO_Grid2DQuad.h"     /* Include 2D quadrilateral block grid header file */
#include "../AMR/QuadTree.h"           /* Include quadtree header file */
#include "../AMR/AMR.h"                /* Include AMR header file */
#include "AdvectDiffuse2DInput.h"      /* Include 2D advection diffusion equation header file */
#include "../ICEM/ICEMCFD.h"           /* Include ICEMCFD header file. */
#include "../System/System_Linux.h"    /* Include System Linux header file. */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */
#include "../HighOrderReconstruction/HighOrder2D.h" /* Include 2D high-order template class header file. */

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D    3

/*!
 * Class: AdvectDiffuse2D_Quad_Block
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
class AdvectDiffuse2D_Quad_Block{
public:

  //! @name Defined public types:
  //@{
  typedef AdvectDiffuse2D_ExactSolutions Exact_Solution_Type;
  typedef AccuracyAssessment2D<AdvectDiffuse2D_Quad_Block> Accuracy_Assessment_Type;
  typedef HighOrder2D<AdvectDiffuse2D_State> HighOrderType; //!< high-order variable data type
  typedef AdvectDiffuse2D_State Soln_State;
  //@}


  //! @name Solution state arrays:
  //@{
  AdvectDiffuse2D_State    **U; //!< Solution state.

  /*! Storage for the solution states at the grid nodes in the calculation of diffusive fluxes */
  AdvectDiffuse2D_State    **U_Nodes;
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

  Grid2D_Quad_Block_HO    Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //! @name Residual and time-stepping arrays:
  //@{
  double                         **dt; //!< Local time step.
  AdvectDiffuse2D_State       ***dUdt; //!< Solution residual.
  AdvectDiffuse2D_State          **Uo; //!< Initial solution.
  static int residual_variable;        //!< Static integer that indicates which variable is used for residual calculations.
  static int Number_of_Residual_Norms; //!< How many Residual norms to plot?
  //@}

  //! @name Solution gradient arrays:
  //@{
  AdvectDiffuse2D_State    **dUdx; //!< Unlimited solution gradient (x-direction).
  AdvectDiffuse2D_State    **dUdy; //!< Unlimited solution gradient (y-direction).
  AdvectDiffuse2D_State     **phi; //!< Solution slope limiter.
  //@}

  //! @name Boundary solution flux arrays:
  //@{
  AdvectDiffuse2D_State   *FluxN, //!< North boundary solution flux.
                          *FluxS, //!< South boundary solution flux.
                          *FluxE, //!< East boundary solution flux.
                          *FluxW; //!< West boundary solution flux.
  //@}

  //! @name Problem indicator flags:
  //@{
  static int      Axisymmetric; //!< Axisymmetric geometry indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static int         Flow_Type; //!< Flow type flag ( required by AMR ).
  static int    Include_Source_Term; //!< Flag for including/excluding the source term in the model equation.
  static int    Include_Advection_Term; //!< Flag for including/excluding the advection term in the model equation.
  static int    Include_Diffusion_Term; //!< Flag for including/excluding the diffusion term in the model equation.
  //@}

  //! @name Boundary condtion reference states:
  //@{
  AdvectDiffuse2D_State   *UoN, //!< Boundary condition reference states for north boundary.
                          *UoS, //!< Boundary condition reference states for south boundary.
                          *UoE, //!< Boundary condition reference states for east boundary.
                          *UoW; //!< Boundary condition reference states for west boundary.
  //! Reference values for north and south boundary conditon reference states
  AdvectDiffuse2D_State Ref_State_BC_North, Ref_State_BC_South; 
  //! Reference values for east and west boundary conditon reference states
  AdvectDiffuse2D_State Ref_State_BC_East, Ref_State_BC_West;
  //@}

  //! @name Accuracy assessment data:
  //@{
  static Exact_Solution_Type *ExactSoln;          //!< Pointer to the exact solution
  static Exact_Solution_Type *ExactSolution(void){ return ExactSoln;} //!< Return the exact solution pointer
  Accuracy_Assessment_Type AssessAccuracy;   //!< Variable to provide access to accuracy assessment subroutines
  //@}
	      
  //! @name Analytically defined boundary conditions:
  //@{
  static AdvectDiffuse2D_InflowField *Inflow;    /*!< Pointer to the inflow field */
  //@}

  //! @name Creation, copy, and assignment constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_Quad_Block(void);

  /* Destructor. */
  ~AdvectDiffuse2D_Quad_Block(void){ deallocate(); }
  //@}

  //! @name Allocate and deallocate functions.
  //@{
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int &Ni, const int &Nj, const int &Ng);

  //! Allocate memory for high-order variables
  void allocate_HighOrder(const int & NumberOfReconstructions,
			  const vector<int> & ReconstructionOrders);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);

  //! Deallocate high-order variable memory for structured quadrilateral solution block.
  void deallocate_HighOrder(void);

  //! Allocate memory for the static memory pool U_Nodes
  void allocate_U_Nodes(const int &_NNi, const int &_NNj);

  //! Deallocate the static memory pool U_Nodes
  void deallocate_U_Nodes(void);
  //@}

  //! Make an identical copy of SolnBlk
  void makeCopy(const AdvectDiffuse2D_Quad_Block &SolnBlk){ *this = SolnBlk; }

  //! @name Bilinear interplation (Zingg & Yarrow).
  //@{
  //! Return solution state at specified node.
  AdvectDiffuse2D_State Un(const int &ii, const int &jj);

  AdvectDiffuse2D_State UnNW(const int &ii, const int &jj); //!< Return solution state at cell NW node.
  AdvectDiffuse2D_State UnNE(const int &ii, const int &jj); //!< Return solution state at cell NE node.
  AdvectDiffuse2D_State UnSE(const int &ii, const int &jj); //!< Return solution state at cell SE node.
  AdvectDiffuse2D_State UnSW(const int &ii, const int &jj); //!< Return solution state at cell SW node.

  //! Return solution state at specified node.
  double un(const int &ii, const int &jj);

  double unNW(const int &ii, const int &jj); //!< Return solution state at cell NW node.
  double unNE(const int &ii, const int &jj); //!< Return solution state at cell NE node.
  double unSE(const int &ii, const int &jj); //!< Return solution state at cell SE node.
  double unSW(const int &ii, const int &jj); //!< Return solution state at cell SW node.
  //@}

  //! @name Field access
  //@{
  const AdvectDiffuse2D_State& U_Node(const int &ii, const int &jj)const { return U_Nodes[ii][jj]; }
  const AdvectDiffuse2D_State& CellSolution(const int &ii, const int &jj)const { return U[ii][jj]; }

  //! @name High-order variables
  //@{
  //! Return the array of high-order variable of the current block
  const HighOrderType * HighOrderVariables(void) const { return HO_Ptr; }
  //! Return the high-order variable in the "Pos" position of the current block
  HighOrderType & HighOrderVariable(const unsigned short int & Pos) { return HO_Ptr[Pos]; }
  const HighOrderType & HighOrderVariable(const unsigned short int & Pos) const { return HO_Ptr[Pos]; }
  //@}

  //@} //end(Field access)

  //! @name Normalization related functions:
  //@{
  //! Get normalization state
  const AdvectDiffuse2D_State getNormalizationState(const int &ii, const int &jj) const { return RefU; }
  //! Set normalization state
  static void Set_Normalization_Reference_State(const AdvectDiffuse2D_State & State){ RefU = State; }
  //@}

  //! @name Evaluate diffusive fluxes.
  //@{
  //! Calculate the solution at mesh nodes
  void Calculate_Nodal_Solutions(void);

  //! Calculate the solution state at a location on an interior inter-cellular face
  void EllipticFluxStateAtInteriorInterface(const int & ii_L, const int & jj_L,
					    const int & ii_R, const int & jj_R,
					    const Vector2D & GQPoint,
					    AdvectDiffuse2D_State & State);

  //! Calculate the solution gradient at the specified interface
  Vector2D InterfaceSolutionGradient(const int & ii_L, const int & jj_L,
				     const int & ii_R, const int & jj_R,
				     const int &Gradient_Reconstruction_Type,
				     const int &Stencil_Flag);

  //! Calculate the solution gradient based on a diamond path reconstruction
  Vector2D DiamondPathGradientReconstruction(const Vector2D &Xl, const AdvectDiffuse2D_State &Ul,
					     const Vector2D &Xd, const AdvectDiffuse2D_State &Ud,
					     const Vector2D &Xr, const AdvectDiffuse2D_State &Ur,
					     const Vector2D &Xu, const AdvectDiffuse2D_State &Uu,
					     const int &stencil_flag);
  //@}

  //! @name Source term member functions
  //@{
  AdvectDiffuse2D_State SourceTerm(const int & ii, const int & jj) const;
  AdvectDiffuse2D_State AxisymmetricSourceTerm(const int & ii, const int & jj) const;
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
  void SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
					       int i_index[], int j_index[],
					       int & n_ptr);
  AdvectDiffuse2D_State PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
							const double &DeltaXToCentroid,
							const double &DeltaYToCentroid) const;
  double PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
					 const double &DeltaXToCentroid,
					 const double &DeltaYToCentroid,
					 const unsigned int &parameter) const;
  AdvectDiffuse2D_State PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
							  const Vector2D &CalculationPoint) const;
  double PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
					   const Vector2D &CalculationPoint,
					   const unsigned int &parameter) const;

  AdvectDiffuse2D_State UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
								 const double &DeltaXToCentroid,
								 const double &DeltaYToCentroid) const;
  AdvectDiffuse2D_State UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
								   const Vector2D &CalculationPoint) const;
  //@}

  //! @name Member functions to compute the inviscid and elliptic flux states at a boundary interface
  //@{
  void InviscidAndEllipticFluxStates_AtBoundaryInterface(const int &BOUNDARY,
							 const int &ii, const int &jj,
							 AdvectDiffuse2D_State &Ul,
							 AdvectDiffuse2D_State &Ur,
							 AdvectDiffuse2D_State &U_face,
							 Vector2D &GradU_face,
							 const int &Gradient_Reconstruction_Type);
  void ViscousFluxStates_AtBoundaryInterface_HighOrder(const int &BOUNDARY,
						       const int &ii, const int &jj,
						       const AdvectDiffuse2D_State &Ul,
						       const AdvectDiffuse2D_State &Ur,
						       AdvectDiffuse2D_State &U_face,
						       const Vector2D &GradUl,
						       const Vector2D &GradUr,
						       Vector2D &GradU_face,
						       const Vector2D &CalculationPoint,
						       const Vector2D &NormalDirection) const;
  void InviscidFluxStates_AtBoundaryInterface_HighOrder(const int &BOUNDARY,
							const int &ii, const int &jj,
							AdvectDiffuse2D_State &Ul,
							AdvectDiffuse2D_State &Ur,
							const Vector2D &CalculationPoint,
							const Vector2D &NormalDirection,
							const unsigned short int Pos = 0) const;
  //@}


  //! @name Member functions to set boundary states
  //@{
  //! Set reference values for boundary reference states
  void Set_Reference_Values_For_Boundary_States(const AdvectDiffuse2D_State & Ref_North,
						const AdvectDiffuse2D_State & Ref_South,
						const AdvectDiffuse2D_State & Ref_East,
						const AdvectDiffuse2D_State & Ref_West);
  //! Set boundary reference states
  void Set_Boundary_Reference_States(void);
  //! Set boundary reference states to default
  void Set_Default_Boundary_Reference_States(void);
  //! Set boundary reference states based on user's input data
  void Set_Boundary_Reference_States_Based_On_Input(const AdvectDiffuse2D_Input_Parameters &IP);
  //@}

  //! @name Residual evaluation functions:
  //@{
  int dUdt_Residual_HighOrder(const AdvectDiffuse2D_Input_Parameters &IP,
			      const int & k_residual,
			      const unsigned short int Pos = 0);
  int dUdt_Residual_Evaluation(const AdvectDiffuse2D_Input_Parameters &IP);
  int dUdt_Residual_Evaluation_HighOrder(const AdvectDiffuse2D_Input_Parameters &IP,
					 const unsigned short int Pos = 0);
  int dUdt_Multistage_Explicit(const int &i_stage,
			       const AdvectDiffuse2D_Input_Parameters &IP);
  int dUdt_Multistage_Explicit_HighOrder(const int &i_stage,
					 const AdvectDiffuse2D_Input_Parameters &IP,
					 const unsigned short int Pos = 0);
  //@}

  //! @name Functions for estimating positivity of elliptic term discretization:
  //@{
  void Analyse_HighOrder_Positivity_For_LaplacianOperator(const int &iCell, const int &jCell,
							  const unsigned short int Pos = 0);
  void Calculate_HighOrder_SolutionCoefficients_For_LaplacianOperator(const int &iCell, const int &jCell,
								      const unsigned short int Pos = 0);
  void Set_HighOrder_InfluenceDomain_For_LaplacianOperator(const int &iCell, const int &jCell,
							   const unsigned short int Pos = 0);
  void Calculate_HighOrder_Discretization_LaplacianOperator(const int &iCell, const int &jCell,
							    const int & k_residual,
							    const unsigned short int Pos = 0);
  void Output_HighOrder_InfluenceDomain_For_LaplacianOperator(ostream &os);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const AdvectDiffuse2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Quad_Block &Soln);
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


  //! @name Member functions used for plotting.
  //@{
  void Output_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				const double &Time,
				const int &Block_Number,
				const int &Output_Title,
				ostream &Out_File,
				const int &IndexHO = 0);

  void Output_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
					   const AdvectDiffuse2D_Input_Parameters &P,
					   const int &Block_Number,
					   const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const AdvectDiffuse2D_Input_Parameters &P,
						 const int &Block_Number,
						 const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const AdvectDiffuse2D_Input_Parameters &P,
						 const int &Block_Number,
						 const int &IndexHO = 0);

  void Output_Tecplot_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
				 const AdvectDiffuse2D_Input_Parameters &P,
				 const int &Block_Number);

  void Output_Cells_Tecplot_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
				       const AdvectDiffuse2D_Input_Parameters &IP,
				       const int &Block_Number);
  //@}

private:
  AdvectDiffuse2D_Quad_Block(const AdvectDiffuse2D_Quad_Block &Soln); //!< Private copy constructor
  AdvectDiffuse2D_Quad_Block & operator = (const AdvectDiffuse2D_Quad_Block &Soln);   //!< Private assignment operator

  int NNi, NNj;		//!< number of nodes in i-direction and j-direction.

  static AdvectDiffuse2D_State RefU;	//!< reference state for normalizing the solution

  //! @name High-order variables and member functions
  //@{
  HighOrderType* HO_Ptr;  //!< pointer to an array of high-order variables
  unsigned short int NumberOfHighOrderVariables; //!< counter for the total number of high-order variables

  //! Allocate high-order variables array.
  void allocate_HighOrder_Array(const int & NumberOfReconstructions);
  //@}

  //! @name Variables for estimating positivity of elliptic term discretization:
  //@{
  IndexType i_index, j_index;	     //!< Storage for indexes of cells that influence the discretization of the Laplace operator.
  std::vector<double>  Laplacian_Coeffs; /*!< The coefficient for each cell in the stencil that
					   appears in the discretization of the Laplace operator. */
  //@}
};


/*!
 * \class SourceTermFunctionalWithPiecewiseLinear
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
					  const int &_jj_,
					  const AdvectDiffuse2D_Quad_Block *_SolnBlk_): ii(_ii_),
											jj(_jj_),
											SolnBlk(_SolnBlk_){ };
  
  //! Return the pointwise value of the source term at a specified location (x,y)
  double operator()(const double &x, const double &y){
    return SolnBlk->U[ii][jj].s(x,y,SolnBlk->PiecewiseLinearSolutionAtLocation(ii,jj,Vector2D(x,y),1));
  }
  
private:
  int ii,jj;	  //!< cell indexes
  const AdvectDiffuse2D_Quad_Block * SolnBlk; //!< pointer to the solution block
  
  //! private default cstr
  SourceTermFunctionalWithPiecewiseLinear();
  //! private assignment operator
  SourceTermFunctionalWithPiecewiseLinear operator=(const SourceTermFunctionalWithPiecewiseLinear &);
};


/*!
 * \class SourceTermFunctionalWithHighOrderPiecewiseInterpolant
 *
 * @brief Definition of the source term functional
 *        with high-order piecewise reconstruction.
 *        The order of this interpolant is runtime dependent.
 *
 * An object of this class calculates the non-linear source 
 * term at a location specified by the Cartesian coordinates (x,y)
 * using information from the cell (ii,jj).
 * The high-order reconstruction of the cell (ii,jj) provides
 * the value of the solution at a given location. Thus, the 
 * source term function can be integrated over the cell area.
 */
class SourceTermFunctionalWithHighOrderPiecewiseInterpolant{
public:
  /*! Constructor with the cell indexes
   *  \param _ii_      the i-index of the cell that provides the information
   *  \param _jj_      the j-index of the cell that provides the information
   *  \param _SolnBlk_ the pointer to the solution block
   *  \param _HO_      the pointer to the high-order object that provides the reconstruction
   */
  SourceTermFunctionalWithHighOrderPiecewiseInterpolant(const int &_ii_,
							const int &_jj_,
							const AdvectDiffuse2D_Quad_Block *_SolnBlk_,
							const AdvectDiffuse2D_Quad_Block::HighOrderType *_HO_): ii(_ii_),
														jj(_jj_),
														SolnBlk(_SolnBlk_),
														HO(_HO_){ };
  
  //! Return the pointwise value of the source term at a specified location (x,y)
  double operator()(const double &x, const double &y){
    return SolnBlk->U[ii][jj].s(x,y,HO->SolutionAtLocation(ii,jj,Vector2D(x,y),1));
  }
  
private:
  int ii,jj;	  //!< cell indexes
  const AdvectDiffuse2D_Quad_Block * SolnBlk; //!< pointer to the solution block
  const AdvectDiffuse2D_Quad_Block::HighOrderType * HO;   //!< pointer to the high-order object
  
  //! private default cstr
  SourceTermFunctionalWithHighOrderPiecewiseInterpolant();
  //! private assignment operator
  SourceTermFunctionalWithHighOrderPiecewiseInterpolant operator=(const SourceTermFunctionalWithHighOrderPiecewiseInterpolant &);
};


/**************************************************************************
 * AdvectDiffuse2D_Quad_Block::Un?? -- Get cell node solution states.     *
 **************************************************************************/
inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block::un?? -- Get cell node solution values.     *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block::unNW(const int &ii, const int &jj) {
  return (un(ii, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block::unNE(const int &ii, const int &jj) {
  return (un(ii+1, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block::unSE(const int &ii, const int &jj) {
  return (un(ii+1, jj));
}

inline double AdvectDiffuse2D_Quad_Block::unSW(const int &ii, const int &jj) {
  return (un(ii, jj));
}

/*********************//**
 * Node solution value. 
 ***********************/
inline double AdvectDiffuse2D_Quad_Block::un(const int &ii, const int &jj) {
  return Un(ii,jj)[1];
}

/***************************************************//**
 * Calculate velocity at the PointOfInterest location
 ****************************************************/
inline Vector2D AdvectDiffuse2D_Quad_Block::VelocityAtLocation(const Vector2D &PointOfInterest) const {
  return (*U)->V(PointOfInterest.x,PointOfInterest.y);
}

/***************************************************//**
 * Calculate velocity at the centroid of cell (ii,jj).
 ****************************************************/
inline Vector2D AdvectDiffuse2D_Quad_Block::VelocityAtCellCentroid(const int & ii, const int & jj) const {
  return U[ii][jj].V(Grid.XCellCentroid(ii,jj), Grid.YCellCentroid(ii,jj));
}

/**************************************************************//**
 * Calculate diffusion coefficient at the centroid of cell (ii,jj).
 * Use the solution stored in the state class as
 * parameter for the diffusion coefficient field.
 ******************************************************************/
inline double AdvectDiffuse2D_Quad_Block::DiffusionCoeffAtCellCentroid(const int & ii, const int & jj) const {
  return U[ii][jj].k(Grid.XCellCentroid(ii,jj), Grid.YCellCentroid(ii,jj), U[ii][jj].u);
}

/**************************************************************//**
 * Determine the stability limit imposed by the source term
 * Use the solution stored in the state class as the parameter 
 * at the centroid of cell (ii,jj) for the non-linear source term field
 ******************************************************************/
inline double AdvectDiffuse2D_Quad_Block::SourceTermStabilityLimit(const int & ii, const int & jj) const {
  return U[ii][jj].SourceTerm->getStabilityLimit(Grid.XCellCentroid(ii,jj),
						 Grid.YCellCentroid(ii,jj),
						 U[ii][jj].u);
}

/**************************************************************//**
 * Determine the axisymmetric source term corresponding to cell (ii,jj).
 * This routine returns actually the integral of the axisymmetric source
 * term over the domain of cell (ii,jj) divided by the area.
 ******************************************************************/
inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block::AxisymmetricSourceTerm(const int & ii,
											const int & jj) const{
  return AdvectDiffuse2D_State(U[ii][jj].s_axi(Grid.Cell[ii][jj].Xc));
}


/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block::evaluate_limiters -- Set flag to evaluate limiters.
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.    *
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block::NumVar -- Returns number of state variables.    *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block::NumVar(void) {
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
inline AdvectDiffuse2D_State 
AdvectDiffuse2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
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
inline double AdvectDiffuse2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
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
inline AdvectDiffuse2D_State 
AdvectDiffuse2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
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
inline double AdvectDiffuse2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
									    const Vector2D &CalculationPoint,
									    const unsigned int &parameter) const{
  return PiecewiseLinearSolutionAtLocation(ii,jj,CalculationPoint)[parameter];
}

/*!
 * Compute the solution provided by the UNLIMITED (i.e. no limiter!) piecewise 
 * linear reconstruction at a location relative to the cell centroid.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param DeltaXToCentroid the X-delta between the point of interest and the centroid
 * \param DeltaYToCentroid the Y-delta between the point of interest and the centroid
 */
inline AdvectDiffuse2D_State 
AdvectDiffuse2D_Quad_Block::UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
								     const double &DeltaXToCentroid,
								     const double &DeltaYToCentroid) const{

  return U[ii][jj] + dUdx[ii][jj]*DeltaXToCentroid + dUdy[ii][jj]*DeltaYToCentroid;
}

/*!
 * Compute the solution provided by the UNLIMITED (i.e. no limiter!) piecewise 
 * linear reconstruction at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline AdvectDiffuse2D_State 
AdvectDiffuse2D_Quad_Block::UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
								       const Vector2D &CalculationPoint) const{

  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	
  return UnlimitedPiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y);
}

/*!
 * Allocate the high-order array.
 */
inline void AdvectDiffuse2D_Quad_Block::allocate_HighOrder_Array(const int & NumberOfReconstructions){

  if (NumberOfReconstructions != NumberOfHighOrderVariables){
    
    // deallocate memory
    deallocate_HighOrder();
    
    // allocate new memory for high-order objects
    HO_Ptr = new HighOrderType [NumberOfReconstructions];
    
    // check for successful memory allocation
    if (HO_Ptr != NULL){
      NumberOfHighOrderVariables = NumberOfReconstructions;
    }
  }
}

/*!
 * Determine the number of neighbouring cells to
 * be used in the piecewise linear reconstruction
 * procedure and set the cell indexes in the 
 * provided containers.
 *
 * \param [in] i the i-index of the reconstructed cell
 * \param [in] j the j-index of the reconstructed cell
 * \param [out] i_index the "i" indexes of the neighbouring cells that form the stencil
 * \param [out] j_index the "j" indexes of the neighbouring cells that form the stencil
 * \param [out] n_ptr the number of neighbouring cells that form the stencil.
 * This number is typically 8, but it can vary for different boundary conditions.
 */
inline void AdvectDiffuse2D_Quad_Block::SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
										int i_index[], int j_index[],
										int &n_pts){

  if (i <= ICl-Nghost || i >= ICu+Nghost || j <= JCl-Nghost || j >= JCu+Nghost) {
    n_pts = 0;
  } else if ((i == ICl-1) && (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
	       Grid.BCtypeW[j] == BC_NEUMANN ||
	       Grid.BCtypeW[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */           
  } else if ((i == ICu+1) && (Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
	       Grid.BCtypeE[j] == BC_NEUMANN ||
	       Grid.BCtypeE[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCl-1) && (Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
	       Grid.BCtypeS[i] == BC_NEUMANN ||
	       Grid.BCtypeS[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCu+1) && (Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
	       Grid.BCtypeN[i] == BC_NEUMANN ||
	       Grid.BCtypeN[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } /* endif */
    } /* endif */
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block External Subroutines.       *
 **************************************************************************/
extern void Write_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk1,
		                const AdvectDiffuse2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Fine,
				  AdvectDiffuse2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block &SolnBlk_Coarse,
				   AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SW,
				   AdvectDiffuse2D_Quad_Block &SolnBlk_Original_SE,
				   AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NW,
				   AdvectDiffuse2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
			   AdvectDiffuse2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
				 AdvectDiffuse2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void ICs(AdvectDiffuse2D_Quad_Block &SolnBlk,
 	        const AdvectDiffuse2D_Input_Parameters &IP,
                AdvectDiffuse2D_State *Wo);

extern void BCs(AdvectDiffuse2D_Quad_Block &SolnBlk,
		const AdvectDiffuse2D_Input_Parameters &IP);

extern double CFL(AdvectDiffuse2D_Quad_Block &SolnBlk,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction(AdvectDiffuse2D_Quad_Block &SolnBlk,
				  const int & Reconstruction_Type,
				  const int & Limiter);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  AdvectDiffuse2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          AdvectDiffuse2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
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

extern int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block &SolnBlk,
				    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters);


/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Multiple Block External Subroutines.     *
 **************************************************************************/

extern AdvectDiffuse2D_Quad_Block* Allocate(AdvectDiffuse2D_Quad_Block *Soln_ptr,
					    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern AdvectDiffuse2D_Quad_Block* Deallocate(AdvectDiffuse2D_Quad_Block *Soln_ptr,
					      AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void ICs(AdvectDiffuse2D_Quad_Block *Soln_ptr,
		AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern void Linear_Reconstruction(AdvectDiffuse2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void BCs(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern double CFL(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr, 
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l1_norm);

extern void L2_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l2_norm);

extern void Max_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      double *max_norm);

extern void Evaluate_Limiters(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

/********************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Multiple Block External Subroutines for Mesh.  *
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

extern int AdvectDiffuse2DQuadSolver(char *Input_File_Name_ptr,
				     int batch_flag);

#endif /* _ADVECTDIFFUSE2D_QUAD_INCLUDED  */

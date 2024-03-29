/*! NavierStokes2DQuad.h: Header file defining 2D Navier-Stokes
                          quadrilateral mesh solution classes. */

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#define _NAVIERSTOKES2D_QUAD_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DState.h"     // Include 2D Navier-Stokes state header file.
#include "NavierStokes2DInput.h"     // Include 2D Navier-Stokes input header file.
#include "../Grid/Cell2D.h"          // Include the 2D cell header file.
#include "../Grid/Grid2DQuad.h"      // Include the 2D quadrilateral multiblock grid header file.
#include "../Grid/HO_Grid2DQuad.h"   /* Include 2D quadrilateral block grid header file */
#include "../AMR/QuadTree.h"         // Include the quadtree header file.
#include "../AMR/AMR.h"              // Include the AMR header file.
#include "../Math/LinearSystems.h"   // Include the linear systems header file.
#include "../ICEM/ICEMCFD.h"         // Include ICEMCFD input header file.
#include "../Turbulent2D/Turbulent2DWallData.h" // Include the turbulent wall data header file.
#include "../System/System_Linux.h"    /* Include System Linux header file. */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */
#include "../HighOrderReconstruction/HighOrder2D.h" /* Include 2D high-order template class header file. */
#include "NavierStokes2D_Cauchy_BCs.h"      /* Include 2D high-order boundary conditions header file,
					       including NavierStokes2D specializations. */
#include "NavierStokes2DExactSolutions.h"   /* Include 2D Navier-Stokes exact solutions header file */

// Define the structures and classes.

#define	NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D  3

// (Parallel) Debugging option.

#define _NS_PARALLEL_DEBUG_

/*!
 * Class: NavierStokes2D_Quad_Block
 *
 * @brief Class definition of the 2D Navier-Stokes solution blocks.
 *
 * \verbatim
 * Member functions
 *       W      -- Return primitive variable solution for the block
 *                 (cell average).
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Wall      -- Return turbulent wall data.
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the solution block residuals.
 *    dWdx      -- Return the unlimited primitive variable solution
 *                 gradients (x-direction) for the block.
 *    dWdy      -- Return the unlimited primitive variable solution
 *                 gradients (y-direction) for the block.
 *    d_dWdx_dW -- Return the gradients of the unlimited primitive 
 *                 variable solution gradients (x-direction) with 
 *                 respect to the primitive variables for the block.
 *    d_dWdy_dW -- Return the gradients of the unlimited primitive 
 *                 variable solution gradients (y-direction) with 
 *                 respect to the primitive variables for the block.
 *     phi      -- Return the solution slope limiters.
 *      Uo      -- Return initial solution state.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in the i-direction 
 *                 (zeta-direction).
 *     ICl      -- Return lower index for cells in the i-direction 
 *                 (zeta-direction).
 *     ICu      -- Return upper index for cells in the i-direction 
 *                 (zeta-direction).
 *     NCj      -- Return number of cells in the j-direction 
 *                 (eta-direction).
 *     JCl      -- Return lower index for cells in the j-direction 
 *                 (eta-direction).
 *     JCu      -- Return upper index for cells in the j-direction 
 *                 (eta-direction).
 *  Nghost      -- Number of ghost cells.
 * Axisymmetric -- Return axisymmetric flow indicator (= 1 for 
 *                 axisymmetric flow, = 0 for planar flow).
 * Compressibility_Effect -- Returns Compressibility corrections indicator
 * 
 * Transition_Model -- Returns the Transition model to be used.
 *     
 * Variable_Prandtl -- Returns if Variable Prandtl number model is to be used
 *                   (=1 for Variable Prandtl on,
 *                   (=0 for Variable Prandtl off).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                   (=1 for limiter freezing on,
 *                   (=0 for limiter freezing off).
 *     WoN      -- Return array of reference states for the
 *                 application of boundary conditions at the north
 *                 boundary of the solution block.
 *     WoS      -- Return array of reference states for the
 *                 application of boundary conditions at the south
 *                 boundary of the solution block.
 *     WoE      -- Return array of reference states for the
 *                 application of boundary conditions at the east
 *                 boundary of the solution block.
 *     WoW      -- Return array of reference states for the
 *                 application of boundary conditions at the west
 *                 boundary of the solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Wn      -- Return primitive variable solution at the
 *                 specified node.
 *      Un      -- Return conserved variable solution at the
 *                 specified node.
 *    WnNW      -- Return primitive variable solution at the
 *                 north-west node.
 *    WnNE      -- Return primitive variable solution at the
 *                 north-east node.
 *    WnSW      -- Return primitive variable solution at the
 *                 south-west node.
 *    WnSE      -- Return primitive variable solution at the
 *                 south-east node.
 *    UnNW      -- Return conserved variable solution at the
 *                 north-west node.
 *    UnNE      -- Return conserved variable solution at the
 *                 north-east node.
 *    UnSW      -- Return conserved variable solution at the
 *                 south-west node.
 *    UnSE      -- Return conserved variable solution at the
 *                 south-east node.
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
 * SubcellReconstruction -- Performs linear subcell solution
 *                          reconstruction used in adaptive mesh
 *                          refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                            conservative flux corrections from fine
 *                            to coarse solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                              conservative flux corrections from
 *                              fine to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D NavierStokes solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class NavierStokes2D_Quad_Block{
public:

  //! @name Defined public types:
  //@{
  typedef NavierStokes2D_ExactSolutions Exact_Solution_Type;
  typedef AccuracyAssessment2D<NavierStokes2D_Quad_Block> Accuracy_Assessment_Type;
  typedef NavierStokes2D_pState Soln_State;
  typedef HighOrder2D<Soln_State> HighOrderType; //!< high-order variable data type. Use primitive variables for reconstruction.
  typedef Cauchy_BCs<Soln_State> BC_Type;
  typedef std::vector<double> DoubleArrayType;
  typedef Grid2D_Quad_Block_HO GridType;
  //@}

  //@{ @name Solution state arrays:
  NavierStokes2D_pState     **W; //!< Primitive solution state.
  NavierStokes2D_cState     **U; //!< Conserved solution state.
  //@}

  //@{ @name Grid block information:
  int                       NCi, //!< Total number of i-direction cells.
                            ICl, //!< First i-direction non-ghost cell counter.
                            ICu; //!< Final i-direction non-ghost cell counter.
  int                       NCj, //!< Total number of j-direction cells.
                            JCl, //!< First j-direction non-ghost cell counter.
                            JCu; //!< Final j-direction non-ghost cell counter.
  int                    Nghost; //!< Number of ghost cells.
  GridType                 Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                   **dt; //!< Local time step.
  NavierStokes2D_cState ***dUdt; //!< Solution residual.
  NavierStokes2D_cState    **Uo; //!< Initial solution state.
  static int  residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  static int Number_of_Residual_Norms;   //!< (default 4 )
  //@}

  //@{ @name Solution gradient arrays:
  NavierStokes2D_pState  **dWdx; //!< Unlimited solution gradient (x-direction).
  NavierStokes2D_pState  **dWdy; //!< Unlimited solution gradient (y-direction).
  NavierStokes2D_pState   **phi; //!< Solution slope limiter.
  //@}

  //@{ @name gradient of Solution gradient w.r.t primitive variables state arrays:
  double                  ***d_dWdx_dW; //!< Unlimited solution gradient (x-direction).
  double                  ***d_dWdy_dW; //!< Unlimited solution gradient (y-direction).
  //@}

  //@{ @name Solution face gradients arrays. For use in viscous NKS preconditioner. Not used in explicit time stepping.
  bool face_grad_arrays_allocated;
  NavierStokes2D_pState       **dWdx_faceN;   //!< north cell face(x-direction).
  NavierStokes2D_pState       **dWdy_faceN;   //!< north cell face(y-direction).
  NavierStokes2D_pState       **dWdx_faceE;   //!< east  cell face(x-direction).
  NavierStokes2D_pState       **dWdy_faceE;   //!< east  cell face(y-direction).
  //@}

  //@{ @name Boundary solution flux arrays:
  NavierStokes2D_cState  *FluxN, //!< North boundary solution flux.
                         *FluxS, //!< South boundary solution flux.
                         *FluxE, //!< East boundary solution flux.
                         *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int                    Axisymmetric; //!< Axisymmetric flow indicator.
  int          Compressibility_Effect; //!< Compressibility Effect indicator
  int                Transition_Model; //!< Transition Model indicator
  int                       Flow_Type; //!< Flow-type flag (inviscid, laminar, or k-omega).
  int                  Freeze_Limiter; //!< Limiter freezing indicator.
  int                Variable_Prandtl; //!< Variable Prandtl number indicator. 
  //@}

  //@{ @name Boundary condition reference states:
  NavierStokes2D_pState    *WoN, //!< Boundary condition reference states for north boundary.
                           *WoS, //!< Boundary condition reference states for south boundary.
                           *WoE, //!< Boundary condition reference states for east boundary.
                           *WoW; //!< Boundary condition reference states for west boundary.
  //! Reference values for north and south boundary conditon reference states
  NavierStokes2D_pState Ref_State_BC_North, Ref_State_BC_South; 
  //! Reference values for east and west boundary conditon reference states
  NavierStokes2D_pState Ref_State_BC_East, Ref_State_BC_West;
  //@}

  //@{ @name Turbulence wall data arrays:
  Turbulent2DWallData    **Wall; //!< Turbulent wall data.
  //@}

  //@{ @name Flow constants:
  Vector2D                Vwall; //!< Specified wall velocity.
  double                  Twall; //!< Specified wall temperature.
  //@}

  //! @name Accuracy assessment data:
  //@{
  static Exact_Solution_Type *ExactSoln;          //!< Pointer to the exact solution
  static Exact_Solution_Type *ExactSolution(void){ return ExactSoln;} //!< Return the exact solution pointer
  Accuracy_Assessment_Type AssessAccuracy;   //!< Variable to provide access to accuracy assessment subroutines
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  NavierStokes2D_Quad_Block(void);

  // Destructor.
  ~NavierStokes2D_Quad_Block(void){ deallocate(); }
  //@}

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Allocate the face gradient arrays. These arrays are only necessary for the implicit (NKS) code.
  void allocate_face_grad_arrays(void);

  //! Allocate memory for high-order variables
  void allocate_HighOrder(const int & NumberOfReconstructions,
			  const vector<int> & ReconstructionOrders,
			  const bool _complete_initialization_ = true);

  //! Allocate memory for high-order boundary conditions
  void allocate_HighOrder_BoundaryConditions(void);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //! Deallocate high-order variable memory for structured quadrilateral solution block.
  void deallocate_HighOrder(void);
  //! Deallocate high-order boundary conditions memory.
  void deallocate_HighOrder_BoundaryConditions(void);
  //@}

  //! Make an identical copy of SolnBlk
  void makeCopy(const NavierStokes2D_Quad_Block &SolnBlk){ *this = SolnBlk; }
  //! @brief Make a copy of the high-order objects
  void copy_HighOrder_Objects(const NavierStokes2D_Quad_Block &SolnBlk);


  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at specified node.
  NavierStokes2D_pState Wn(const int ii, const int jj);

  //! Return conserverd solution state at specified node.
  NavierStokes2D_cState Un(const int ii, const int jj);

  NavierStokes2D_pState WnNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  NavierStokes2D_pState WnNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  NavierStokes2D_pState WnSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  NavierStokes2D_pState WnSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  NavierStokes2D_cState UnNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  NavierStokes2D_cState UnNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  NavierStokes2D_cState UnSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  NavierStokes2D_cState UnSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.
  //@}

  //@{ @name Bilinear interplation (Connell & Holmes).
  //! Return primitive solution state at specified node.
  NavierStokes2D_pState Ww(const int ii, const int jj);

  //! Retern conserverd solution state at specified node.
  NavierStokes2D_cState Uw(const int ii, const int jj);

  //! Return primitive solution state at cell nodes.
  NavierStokes2D_pState WwNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  NavierStokes2D_pState WwNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  NavierStokes2D_pState WwSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  NavierStokes2D_pState WwSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  //! Return conserved solution state at cell nodes.
  NavierStokes2D_cState UwNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  NavierStokes2D_cState UwNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  NavierStokes2D_cState UwSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  NavierStokes2D_cState UwSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.
  //@}

  //! @name Field access
  //@{
  const NavierStokes2D_pState& CellSolution(const int &ii, const int &jj) const { return W[ii][jj]; }

  //! @name High-order variables
  //@{
  //! Return the array of high-order variable of the current block
  const HighOrderType * HighOrderVariables(void) const { return HO_Ptr; }
  //! Return the high-order variable in the "Pos" position of the current block
  HighOrderType & HighOrderVariable(const unsigned short int & Pos) { return HO_Ptr[Pos]; }
  const HighOrderType & HighOrderVariable(const unsigned short int & Pos) const { return HO_Ptr[Pos]; }
  //! Return the number of high-order variables
  const unsigned short int & NumberOfHighOrderObjects(void) const { return NumberOfHighOrderVariables; }
  //@}

  //! @name High-order boundary conditions (used mostly for constrained reconstruction)
  //@{
  const BC_Type * BC_NorthCell(void) { return HO_WoN;}
  BC_Type & BC_NorthCell(const int &iCell){ return HO_WoN[iCell]; }
  const BC_Type & BC_NorthCell(const int &iCell) const { return HO_WoN[iCell]; }
  const BC_Type * BC_SouthCell(void) { return HO_WoS;}
  BC_Type & BC_SouthCell(const int &iCell){ return HO_WoS[iCell]; }
  const BC_Type & BC_SouthCell(const int &iCell) const { return HO_WoS[iCell]; }
  const BC_Type * BC_EastCell(void) { return HO_WoE;}
  BC_Type & BC_EastCell(const int &jCell){ return HO_WoE[jCell]; }
  const BC_Type & BC_EastCell(const int &jCell) const { return HO_WoE[jCell]; }
  const BC_Type * BC_WestCell(void) { return HO_WoW;}
  BC_Type & BC_WestCell(const int &jCell){ return HO_WoW[jCell]; }
  const BC_Type & BC_WestCell(const int &jCell) const { return HO_WoW[jCell]; }

  void BCs_HighOrder(void);
  //@}

  //@} //end(Field access)

  //! @name Normalization related functions:
  //@{
  //! Get normalization state
  const NavierStokes2D_pState getNormalizationState(const int &ii, const int &jj) const { return RefW; }
  //! Set the normalization state which is used in the smoothness indicator computation.
  static void Set_Normalization_Reference_State(const NavierStokes2D_pState & State){ 
    RefW = NavierStokes2D_pState(State.rho, State.a(), State.a(), State.p, State.k, State.omega, State.ke, State.ee); 
  }
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //! @name Member functions to compute the piecewise linear solution at a particular location
  //@{
  void SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
					       int i_index[], int j_index[],
					       int & n_ptr);
  NavierStokes2D_pState PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
							const double &DeltaXToCentroid,
							const double &DeltaYToCentroid) const;
  double PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
					 const double &DeltaXToCentroid,
					 const double &DeltaYToCentroid,
					 const unsigned int &parameter) const;
  NavierStokes2D_pState PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
							  const Vector2D &CalculationPoint) const;
  double PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
					   const Vector2D &CalculationPoint,
					   const unsigned int &parameter) const;

  NavierStokes2D_pState UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
							  const double &DeltaXToCentroid,
							  const double &DeltaYToCentroid) const;
  NavierStokes2D_pState UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
							    const Vector2D &CalculationPoint) const;
  //! @brief Member function to compute the solution entropy at a particular location
  double SolutionEntropyAtLocation(const int &ii, const int &jj,
				   const Vector2D &CalculationPoint,
				   const unsigned int &parameter = 0) const;
  //! @brief Member function to compute the solution pressure at a particular location
  double SolutionPressureAtCoordinates(const int &ii, const int &jj,
				       const double &X_Coord, const double &Y_Coord) const;
  //! @brief Member function to compute the wall shear stress at a particular location
  double WallShearStress_HighOrder(const int &ii, const int &jj,
				   const Vector2D &CalculationPoint,
				   const Vector2D &normal_dir) const;
  //@}

  //! @name Member functions to set boundary states
  //@{
  //! Set reference values for boundary reference states
  void Set_Reference_Values_For_Boundary_States(const NavierStokes2D_pState & Ref_North,
						const NavierStokes2D_pState & Ref_South,
						const NavierStokes2D_pState & Ref_East,
						const NavierStokes2D_pState & Ref_West);
  //! Set boundary reference states
  void Set_Boundary_Reference_States(const bool &SetNorth,
				     const bool &SetSouth,
				     const bool &SetEast,
				     const bool &SetWest);
  //! Set boundary reference states to default
  void Set_Default_Boundary_Reference_States(void);
  //! Set boundary reference states based on user's input data
  void Set_Boundary_Reference_States_Based_On_Input(const NavierStokes2D_Input_Parameters &IP);
  //! @brief Set physical boundary condition constraints based on the current flow state and the BC_Type
  void EnsurePhysicalBCsConstraints(const int & BOUNDARY, const int & BndCellIndex);
  //@}

  //! @name Residual evaluation functions:
  //@{
  int dUdt_Residual_Evaluation(NavierStokes2D_Input_Parameters &IP);
  int dUdt_Residual_HighOrder(const NavierStokes2D_Input_Parameters &IP,
			      const int & k_residual,
			      const bool & UseTimeStep,
			      const unsigned short int Pos = 0);
  int dUdt_Residual_Evaluation_HighOrder(const NavierStokes2D_Input_Parameters &IP,
					 const unsigned short int Pos = 0);
  int dUdt_Multistage_Explicit(const int &i_stage,
			       NavierStokes2D_Input_Parameters &IP);
  int dUdt_Multistage_Explicit_HighOrder(const int &i_stage,
					 const NavierStokes2D_Input_Parameters &IP,
					 const unsigned short int Pos = 0);
  //@}

  //! @name Variables/functions for AMR:
  //@{
  //! Return the array of refinement criteria
  DoubleArrayType & Refinement_Criteria(void) { return refinement_criteria; }
  //! Return the array of refinement criteria as constant
  const DoubleArrayType & Refinement_Criteria(void) const { return refinement_criteria; }
  //! Return a particular refinement criterion in the refinement criteria array
  double & Refinement_Criterion(const int & index) {return refinement_criteria[index]; }
  //! Return a particular refinement criterion in the refinement criteria array as constant
  const double & Refinement_Criterion(const int & index) const {return refinement_criteria[index]; }
  //! Return number of refinement criteria
  int Number_Of_Refinement_Criteria(void) const { return refinement_criteria.size(); }
  void Calculate_Refinement_Criteria_HighOrder(double *refinement_criteria,
					       NavierStokes2D_Input_Parameters &IP,
					       int &number_refinement_criteria);
  //@}

  //! @name Functions for flux calculation:
  //@{
  //! @brief Calculate Riemann flux
  NavierStokes2D_cState RiemannFlux_n(const int & Flux_Function,
				      const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr,
				      const Vector2D &normal_dir) const;
  //! @brief Check the validity of the state
  void Validate_Primitive_SolnState(NavierStokes2D_pState & W,
				    const int &iCell,
				    const int &jCell,
				    const std::string &Ref,
				    const int &IndexHO) const;
  void InviscidFluxStates_AtBoundaryInterface_HighOrder(const int &BOUNDARY,
							const int &ii, const int &jj,
							NavierStokes2D_pState &Wl,
							NavierStokes2D_pState &Wr,
							const Vector2D &CalculationPoint,
							const Vector2D &NormalDirection,
							const unsigned short int Pos = 0) const;
  void ViscousFluxStates_AtBoundaryInterface_HighOrder(const int &BOUNDARY,
						       const int &ii, const int &jj,
						       const NavierStokes2D_pState &Wl, const NavierStokes2D_pState &Wr,
						       NavierStokes2D_pState &W_face,
						       const NavierStokes2D_pState &dWdxL, const NavierStokes2D_pState &dWdyL,
						       const NavierStokes2D_pState &dWdxR, const NavierStokes2D_pState &dWdyR,
						       NavierStokes2D_pState &dWdx_face,
						       NavierStokes2D_pState &dWdy_face,
						       const Vector2D &CalculationPoint,
						       const Vector2D &NormalDirection);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const NavierStokes2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, NavierStokes2D_Quad_Block &Soln);

  void Output_Flat_Plate_Tecplot_HighOrder(const int Block_Number,
					   const int Output_Title_Soln,
					   ostream &Out_File_Soln,
					   const int Output_Title_Skin,
					   ostream &Out_File_Skin,
					   const NavierStokes2D_pState &Winf,
					   const double &plate_length,
					   double &l1_norm,
					   double &l2_norm,
					   double &max_norm,
					   double &area,
					   int &numberofactivecells,
					   double &l1_norm_cf,
					   double &l2_norm_cf,
					   double &max_norm_cf,
					   double &area_cf,
					   int &numberofactivecells_cf);
  //@}

  //@{ @name Member functions required for message passing.
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
  //! Linear subcell solution reconstruction within given computational cell.
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

#ifdef _NS_PARALLEL_DEBUG_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  static ofstream dout;

  int Open_Debug_Output_File(const int &ThisCPU) {
    strcpy(output_file_name,"debug");
    strcat(output_file_name,"_cpu");
    sprintf(extension,"%.6d",ThisCPU);
    strcat(extension,".txt");
    strcat(output_file_name,extension);
    output_file_name_ptr = output_file_name;
    dout.open(output_file_name_ptr,ios::out);
    if (dout.bad()) return 1;
    return 0;
  }

  void Close_Debug_Output_File(void) {
    dout.close();
  }
#endif

  //! @name Member functions used for plotting.
  //@{
  void Output_Nodes_Tecplot_HighOrder(const NavierStokes2D_Input_Parameters &IP,
				      const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int & StartI_CellIndex,
				      const int & EndI_CellIndex,
				      const int & StartJ_CellIndex,
				      const int & EndJ_CellIndex,
				      const int &IndexHO = 0);

  void Output_Tecplot_HighOrder(const NavierStokes2D_Input_Parameters &IP,
				const int &Number_of_Time_Steps,
				const double &Time,
				const int &Block_Number,
				const int &Output_Title,
				ostream &Out_File,
				const int &IndexHO = 0);

  void Output_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
					   const NavierStokes2D_Input_Parameters &IP,
					   const int &Block_Number,
					   const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder(const NavierStokes2D_Input_Parameters &IP,
				      const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Nodes_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const NavierStokes2D_Input_Parameters &IP,
						 const int &Block_Number,
						 const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder(const int &Number_of_Time_Steps,
				      const double &Time,
				      const int &Block_Number,
				      const int &Output_Title,
				      ostream &Out_File,
				      const int &IndexHO = 0);

  void Output_Cells_Tecplot_HighOrder_Debug_Mode(AdaptiveBlock2D_List &Soln_Block_List,
						 const NavierStokes2D_Input_Parameters &P,
						 const int &Block_Number,
						 const int &IndexHO = 0);
  //@}

private:
  NavierStokes2D_Quad_Block(const NavierStokes2D_Quad_Block &Soln); //!< Private copy constructor.
  NavierStokes2D_Quad_Block & operator = (const NavierStokes2D_Quad_Block &Soln);   //!< Private assignment operator.
  
  static NavierStokes2D_pState RefW;	//!< reference state for normalizing the solution

  //! @name High-order variables and member functions
  //@{
  HighOrderType* HO_Ptr;  //!< pointer to an array of high-order variables
  unsigned short int NumberOfHighOrderVariables; //!< counter for the total number of high-order variables

  //! Allocate high-order variables array.
  void allocate_HighOrder_Array(const int & NumberOfReconstructions);
  //@}

  //! @name High-order boundary conditions (used mostly for constrained reconstruction)
  //@{
  BC_Type  *HO_WoN, 		//!< High-order boundary condition reference states for North boundary
    *HO_WoS,            	//!< High-order boundary condition reference states for South boundary
    *HO_WoE,			//!< High-order boundary condition reference states for East boundary
    *HO_WoW;            	//!< High-order boundary condition reference states for West boundary
  //@}

  //! @name Variables/functions for AMR:
  //@{
  DoubleArrayType refinement_criteria;
  //@}
};

//! Default constructor.
inline NavierStokes2D_Quad_Block::NavierStokes2D_Quad_Block(void):
  AssessAccuracy(this),
  Ref_State_BC_North(0.0), Ref_State_BC_South(0.0),
  Ref_State_BC_East(0.0), Ref_State_BC_West(0.0),
  HO_Ptr(NULL), NumberOfHighOrderVariables(0)
{
  // Problem flags:
  Axisymmetric   = OFF;
  Compressibility_Effect = OFF;
  Transition_Model = OFF;
  Flow_Type      = FLOWTYPE_INVISCID;
  Freeze_Limiter = OFF;
  Variable_Prandtl = OFF;
  Vwall = Vector2D_ZERO;
  Twall = ZERO;
  // Grid size and variables:
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
  // Solution variables:
  W = NULL; U = NULL;
  dt = NULL; dUdt = NULL; Uo = NULL;
  dWdx = NULL; dWdy = NULL; phi = NULL;
  d_dWdx_dW = NULL; d_dWdy_dW = NULL;
  face_grad_arrays_allocated = false;
  dWdx_faceN = NULL; dWdy_faceN = NULL;
  dWdx_faceE = NULL; dWdy_faceE = NULL;
  FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
  WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
  HO_WoN = NULL; HO_WoS = NULL; HO_WoE = NULL; HO_WoW = NULL;
  // Turbulent wall data:
  Wall = NULL;

  // Get access to the NavierStokes2D_ExactSolutions object
  ExactSoln = &NavierStokes2D_ExactSolutions::getInstance();

}

//! Private copy constructor. (shallow copy)
inline NavierStokes2D_Quad_Block::NavierStokes2D_Quad_Block(const NavierStokes2D_Quad_Block &Soln):
  AssessAccuracy(this), HO_Ptr(NULL)
{
  // Problem flags:
  Axisymmetric   = Soln.Axisymmetric;
  Compressibility_Effect = Soln.Compressibility_Effect;
  Transition_Model = Soln.Transition_Model;
  Flow_Type      = Soln.Flow_Type;
  Freeze_Limiter = Soln.Freeze_Limiter;
  Variable_Prandtl = Soln.Variable_Prandtl;
  Vwall = Soln.Vwall;
  Twall = Soln.Twall;
  // Grid size and variables:
  NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu;
  NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu;
  Nghost = Soln.Nghost;
  Grid  = Soln.Grid;
  // Solution variables:
  W = Soln.W; U = Soln.U;
  dt = Soln.dt; dUdt = Soln.dUdt; Uo = Soln.Uo;
  dWdx = Soln.dWdx; dWdy = Soln.dWdy; phi = Soln.phi;
  d_dWdx_dW = Soln.d_dWdx_dW; d_dWdy_dW = Soln.d_dWdy_dW;
  face_grad_arrays_allocated = Soln.face_grad_arrays_allocated;
  dWdx_faceN = Soln.dWdx_faceN; dWdy_faceN = Soln.dWdy_faceN;
  dWdx_faceE = Soln.dWdx_faceE; dWdy_faceE = Soln.dWdy_faceE;
  FluxN = Soln.FluxN; FluxS = Soln.FluxS; 
  FluxE = Soln.FluxE; FluxW = Soln.FluxW;
  WoN   = Soln.WoN;   WoS   = Soln.WoS;
  WoE   = Soln.WoE;   WoW   = Soln.WoW;
  HO_WoN = Soln.HO_WoN; HO_WoS = Soln.HO_WoS; HO_WoE = Soln.HO_WoE; HO_WoW = Soln.HO_WoW;
  // Turbulent wall data:
  Wall = Soln.Wall;

  Ref_State_BC_North = Soln.Ref_State_BC_North;
  Ref_State_BC_South = Soln.Ref_State_BC_South;
  Ref_State_BC_East = Soln.Ref_State_BC_East;
  Ref_State_BC_West = Soln.Ref_State_BC_West;

  HO_Ptr = Soln.HO_Ptr;
  NumberOfHighOrderVariables = Soln.NumberOfHighOrderVariables;
}


/**********************************************************************
 * NavierStokes2D_Quad_Block::allocate -- Allocate memory.            *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  int i,j,k,m;
  assert(Ni > 1 && Nj > 1 && Ng > 1);

  // Check to see if the current block dimensions differ from the required ones.
  if ( (Nghost != Ng) || (NCi != Ni+2*Ng) || (NCj != Nj+2*Ng) ){ 
    
    // free the memory if there is memory allocated
    deallocate();

    // allocate new memory
    Grid.allocate(Ni,Nj,Ng);
    NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
    NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
    W = new NavierStokes2D_pState*[NCi]; U = new NavierStokes2D_cState*[NCi];
    dt = new double*[NCi]; dUdt = new NavierStokes2D_cState**[NCi];
    dWdx = new NavierStokes2D_pState*[NCi]; dWdy = new NavierStokes2D_pState*[NCi];
    phi = new NavierStokes2D_pState*[NCi]; 
    d_dWdx_dW = new double **[NCi];
    d_dWdy_dW = new double **[NCi];
    Uo = new NavierStokes2D_cState*[NCi];
    Wall = new Turbulent2DWallData*[NCi];
    for (i = 0; i < NCi; i++) {
      W[i] = new NavierStokes2D_pState[NCj]; U[i] = new NavierStokes2D_cState[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new NavierStokes2D_cState*[NCj];
      d_dWdx_dW[i] = new double  *[NCj];
      d_dWdy_dW[i] = new double  *[NCj];
      for ( j = 0; j < NCj; j++){
	dUdt[i][j] = new NavierStokes2D_cState[NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D];
	d_dWdx_dW[i][j] = new double [5];
	d_dWdy_dW[i][j] = new double [5];
      }
      dWdx[i] = new NavierStokes2D_pState[NCj]; dWdy[i] = new NavierStokes2D_pState[NCj];
      phi[i] = new NavierStokes2D_pState[NCj]; Uo[i] = new NavierStokes2D_cState[NCj];
      Wall[i] = new Turbulent2DWallData[NCj];
    } /* endfor */
    WoN = new NavierStokes2D_pState[NCi]; WoS = new NavierStokes2D_pState[NCi];
    WoE = new NavierStokes2D_pState[NCj]; WoW = new NavierStokes2D_pState[NCj];
    FluxN = new NavierStokes2D_cState[NCi]; FluxS = new NavierStokes2D_cState[NCi];
    FluxE = new NavierStokes2D_cState[NCj]; FluxW = new NavierStokes2D_cState[NCj];
    // Set the solution residuals, gradients, limiters, and other values to zero.
    for ( j = JCl-Nghost; j <= JCu+Nghost; j++) {
      for ( i = ICl-Nghost; i <= ICu+Nghost; i++) {
	for ( k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; k++){
	  dUdt[i][j][k].Vacuum();
	} /* endfor */
	dWdx[i][j].Vacuum(); dWdy[i][j].Vacuum();
	phi[i][j].Vacuum(); Uo[i][j].Vacuum();
	dt[i][j] = ZERO;
      } /* endfor */
    } /* endfor */

  }/* endif */
}

inline void NavierStokes2D_Quad_Block::allocate_face_grad_arrays(void) {
  if (face_grad_arrays_allocated) { return; }
  
  // call allocate() then call this
  assert(NCi > 0 && NCj > 0);

  dWdx_faceN = new NavierStokes2D_pState*[NCi]; dWdy_faceN = new NavierStokes2D_pState*[NCi];
  dWdx_faceE = new NavierStokes2D_pState*[NCi]; dWdy_faceE = new NavierStokes2D_pState*[NCi];
  for (int i = 0; i < NCi; i++) {
    dWdx_faceN[i] = new NavierStokes2D_pState[NCj]; dWdy_faceN[i] = new NavierStokes2D_pState[NCj];
    dWdx_faceE[i] = new NavierStokes2D_pState[NCj]; dWdy_faceE[i] = new NavierStokes2D_pState[NCj];

    for (int j = 0; j < NCj; j++) {
      dWdx_faceN[i][j].Vacuum(); dWdy_faceN[i][j].Vacuum();
      dWdx_faceE[i][j].Vacuum(); dWdy_faceE[i][j].Vacuum();
    }
  }

  face_grad_arrays_allocated = true;
}

/*****************************************//**
 * Allocate memory for high-order variables.
 ********************************************/
inline void NavierStokes2D_Quad_Block::allocate_HighOrder(const int & NumberOfReconstructions,
							  const vector<int> & ReconstructionOrders,
							  const bool _complete_initialization_){

  bool _pseudo_inverse_allocation_(false);
  int i;

  // Decide whether to allocate the pseudo-inverse
  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT){
    _pseudo_inverse_allocation_ = true;
  }

  // Re-allocate new memory if necessary
  if (NumberOfReconstructions != NumberOfHighOrderVariables){

    // allocate the high-order array
    allocate_HighOrder_Array(NumberOfReconstructions);
    
    // set the reconstruction order of each high-order object
    for (i=0; i<NumberOfHighOrderVariables; ++i){
      if (_complete_initialization_){
	// initialize the high-order variable completely 
	HO_Ptr[i].InitializeVariable(ReconstructionOrders[i],
				     Grid,
				     _pseudo_inverse_allocation_);
      } else {
	// initialize the basic high-order variable
	HO_Ptr[i].InitializeBasicVariable(ReconstructionOrders[i],
					  Grid,
					  _pseudo_inverse_allocation_);
      }
    }

  } else {
    // check the reconstruction orders
    for (i=0; i<ReconstructionOrders.size(); ++i){
      if (HighOrderVariable(i).RecOrder() != ReconstructionOrders[i]){
	// change the reconstruction order of the high-order object
	HO_Ptr[i].SetReconstructionOrder(ReconstructionOrders[i]);
      }
    } // endfor
  }// endif
}

/***************************************************//**
 * Allocate memory for high-order boundary conditions
 * if the corresponding boundary reconstruction is 
 * constrained.
 * Assume that Grid and block indexes have been setup!
 *******************************************************/
inline void NavierStokes2D_Quad_Block::allocate_HighOrder_BoundaryConditions(void){
  
  int i,j;
  
  // allocate North BCs
  if ( Grid.IsWestExtendNorthBoundaryReconstructionConstrained() ||  
       Grid.IsNorthBoundaryReconstructionConstrained() || 
       Grid.IsEastExtendNorthBoundaryReconstructionConstrained() ){
    
    if (HO_WoN != NULL){
      // deallocate memory
      delete [] HO_WoN; HO_WoN = NULL;
    }

    // allocate new memory
    HO_WoN = new BC_Type[NCi];

    // allocate BC memory for each constrained Gauss quadrature point
    for (i=0; i<NCi; ++i){
      BC_NorthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_North(i,JCu),
					  Grid.BCtypeN[i]);
    }

  } else if ( HO_WoN != NULL){
    // deallocate memory
    delete [] HO_WoN; HO_WoN = NULL;
  }

  // allocate South BCs
  if ( Grid.IsWestExtendSouthBoundaryReconstructionConstrained() ||
       Grid.IsSouthBoundaryReconstructionConstrained() || 
       Grid.IsEastExtendSouthBoundaryReconstructionConstrained() ){

    if (HO_WoS != NULL){
      // deallocate memory
      delete [] HO_WoS; HO_WoS = NULL;
    }

    // allocate new memory    
    HO_WoS = new BC_Type[NCi];

    // allocate BC memory for each constrained Gauss quadrature point
    for (i=0; i<NCi; ++i){
      BC_SouthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_South(i,JCl),
					  Grid.BCtypeS[i]);
    }    
  } else if (HO_WoS != NULL){
    // deallocate memory
    delete [] HO_WoS; HO_WoS = NULL;
  }

  // allocate East BCs
  if ( Grid.IsSouthExtendEastBoundaryReconstructionConstrained() ||
       Grid.IsEastBoundaryReconstructionConstrained() ||
       Grid.IsNorthExtendEastBoundaryReconstructionConstrained() ){

    if (HO_WoE != NULL){
      // deallocate memory
      delete [] HO_WoE; HO_WoE = NULL;
    }

    // allocate new memory    
    HO_WoE = new BC_Type[NCj];

    // allocate BC memory for each constrained Gauss quadrature point
    for (j=0; j<NCj; ++j){
      BC_EastCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_East(ICu,j),
					 Grid.BCtypeE[j]);
    }
  } else if (HO_WoE != NULL){
    // deallocate memory
    delete [] HO_WoE; HO_WoE = NULL;
  }

  // allocate West BCs
  if ( Grid.IsSouthExtendWestBoundaryReconstructionConstrained() ||
       Grid.IsWestBoundaryReconstructionConstrained() || 
       Grid.IsNorthExtendWestBoundaryReconstructionConstrained() ){

    if (HO_WoW != NULL){
      // deallocate memory
      delete [] HO_WoW; HO_WoW = NULL;
    }

    // allocate new memory    
    HO_WoW = new BC_Type[NCj];

    // allocate BC memory for each constrained Gauss quadrature point
    for (j=0; j<NCj; ++j){
      BC_WestCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_West(ICl,j),
					 Grid.BCtypeW[j]);
    }

  } else if (HO_WoW != NULL){
    // deallocate memory
    delete [] HO_WoW; HO_WoW = NULL;
  }

}


/**********************************************************************
 * NavierStokes2D_Quad_Block::deallocate -- Deallocate memory.        *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::deallocate(void) {
  if (U != NULL){
    /* free the memory if there is memory allocated */
    int i,j;
    Grid.deallocate();
    for (i = 0; i < NCi; i++) {
      delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL; 
      for (j = 0; j < NCj; j++) { 
	delete []dUdt[i][j]; dUdt[i][j] = NULL; 
	delete []d_dWdx_dW[i][j]; d_dWdx_dW[i][j] = NULL;
	delete []d_dWdy_dW[i][j]; d_dWdy_dW[i][j] = NULL;
      }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
      delete []d_dWdx_dW[i]; d_dWdx_dW[i] = NULL;
      delete []d_dWdy_dW[i]; d_dWdy_dW[i] = NULL;
      delete []phi[i];  phi[i]  = NULL; delete []Uo[i]; Uo[i] = NULL;
      delete []Wall[i]; Wall[i] = NULL; 
    } /* endfor */
    delete []W; W = NULL; delete []U; U = NULL;
    delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
    delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL; 
    delete []d_dWdx_dW; d_dWdx_dW = NULL;
    delete []d_dWdy_dW; d_dWdy_dW = NULL;
    delete []phi; phi = NULL; delete []Uo; Uo = NULL;
    delete []Wall; Wall = NULL;
    delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
    delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
    delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
    delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
    if (face_grad_arrays_allocated) {
      for (i = 0; i < NCi; i++) {
	delete []dWdx_faceN[i]; delete []dWdy_faceN[i]; 
	delete []dWdx_faceE[i]; delete []dWdy_faceE[i]; 
      }
      delete []dWdx_faceN; delete []dWdy_faceN; 
      delete []dWdx_faceE; delete []dWdy_faceE; 

      face_grad_arrays_allocated = false;
    }
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;

    deallocate_HighOrder();
    deallocate_HighOrder_BoundaryConditions();   
  } // endif
}

/******************************************//**
 * Deallocate memory for high-order variables
 *********************************************/
inline void NavierStokes2D_Quad_Block::deallocate_HighOrder(void) {
  delete []HO_Ptr; HO_Ptr = NULL;
  NumberOfHighOrderVariables = 0;
}

/****************************************************//**
 * Deallocate memory for high-order boundary conditions
 *******************************************************/
inline void NavierStokes2D_Quad_Block::deallocate_HighOrder_BoundaryConditions(void) {
  if (HO_WoN != NULL){
    delete [] HO_WoN; HO_WoN = NULL;
  }
  if (HO_WoS != NULL){
    delete [] HO_WoS; HO_WoS = NULL;
  }
  if (HO_WoE != NULL){
    delete [] HO_WoE; HO_WoE = NULL;
  }
  if (HO_WoW != NULL){
    delete [] HO_WoW; HO_WoW = NULL;
  }
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Wn -- Node primitive variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::Wn(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
         eta1, zeta1, eta2, zeta2, eta, zeta;
  NavierStokes2D_pState A, B, C, D;
  x  = Grid.Node[ii][jj].X.x;
  y  = Grid.Node[ii][jj].X.y;
  ax = Grid.Cell[ii-1][jj-1].Xc.x;
  bx = Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x;
  cx = Grid.Cell[ii  ][jj-1].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x;
  dx = Grid.Cell[ii  ][jj  ].Xc.x + Grid.Cell[ii-1][jj-1].Xc.x -
       Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii  ][jj-1].Xc.x;
  ay = Grid.Cell[ii-1][jj-1].Xc.y;
  by = Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y;
  cy = Grid.Cell[ii  ][jj-1].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y;
  dy = Grid.Cell[ii  ][jj  ].Xc.y + Grid.Cell[ii-1][jj-1].Xc.y -
       Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by+dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) {
      zeta1 = -cc/bb;
    } else {
      zeta1 = -cc/sgn(bb)*(TOLER*TOLER);
    }
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    } else {
      eta1 = HALF;
    }
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) {
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa;
    } else { zeta1 = -HALF*bb/aa;
    }
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    }
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) {
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa;
    } else {
      zeta2 = -HALF*bb/aa;
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) {
      eta2 = -ONE;
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  A = W[ii-1][jj-1];
  B = W[ii-1][jj  ] - W[ii-1][jj-1];
  C = W[ii  ][jj-1] - W[ii-1][jj-1];
  D = W[ii  ][jj  ] + W[ii-1][jj-1] - W[ii-1][jj  ] - W[ii  ][jj-1];
  return A + B*zeta + C*eta + D*zeta*eta;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Un -- Node conserved variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::Un(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  NavierStokes2D_cState A, B, C, D;
  x  = Grid.Node[ii][jj].X.x;
  y  = Grid.Node[ii][jj].X.y;
  ax = Grid.Cell[ii-1][jj-1].Xc.x;
  bx = Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x; 
  cx = Grid.Cell[ii  ][jj-1].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x; 
  dx = Grid.Cell[ii  ][jj  ].Xc.x + Grid.Cell[ii-1][jj-1].Xc.x -
       Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii  ][jj-1].Xc.x;
  ay = Grid.Cell[ii-1][jj-1].Xc.y;
  by = Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y; 
  cy = Grid.Cell[ii  ][jj-1].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y; 
  dy = Grid.Cell[ii  ][jj  ].Xc.y + Grid.Cell[ii-1][jj-1].Xc.y -
       Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by + dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }
    else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
    else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { 
      zeta1 = -HALF*bb/aa; } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else { 
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } 
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa;
    } else {
      zeta2 = -HALF*bb/aa; 
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) {
      eta2 = -ONE; 
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  A = U[ii-1][jj-1];
  B = U[ii-1][jj  ] - U[ii-1][jj-1]; 
  C = U[ii  ][jj-1] - U[ii-1][jj-1];
  D = U[ii  ][jj  ] + U[ii-1][jj-1] - U[ii-1][jj  ] - U[ii  ][jj-1];
  return A+B*zeta+C*eta+D*zeta*eta;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Wn?? -- Cell node primitive solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnNW(const int ii, const int jj) {
  return Wn(ii,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnNE(const int ii, const int jj) {
  return Wn(ii+1,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnSE(const int ii, const int jj) {
  return Wn(ii+1,jj);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnSW(const int ii, const int jj) {
  return Wn(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Un?? -- Cell node conserved solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnNW(const int ii, const int jj) {
  return Un(ii,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnNE(const int ii, const int jj) {
  return Un(ii+1,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnSE(const int ii, const int jj) {
  return Un(ii+1,jj);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnSW(const int ii, const int jj) {
  return Un(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Ww -- Node primitive variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::Ww(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  NavierStokes2D_pState W1, W2, W3, W4;
  // Summarize cell-centres and states.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc; W1 = W[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc; W2 = W[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc; W3 = W[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc; W4 = W[ii  ][jj  ];
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) + sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) + sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = ONE + lambda*(X1 - X0);
  w2 = ONE + lambda*(X2 - X0);
  w3 = ONE + lambda*(X3 - X0);
  w4 = ONE + lambda*(X4 - X0);
  // Determine the interpolated state:
  return (w1*W1 + w2*W2 + w3*W3+ w4*W4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Uw -- Node conserved variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::Uw(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  NavierStokes2D_cState U1, U2, U3, U4;
  // Summarize cell-centres and state.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc; U1 = U[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc; U2 = U[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc; U3 = U[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc; U4 = U[ii  ][jj  ];
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) + sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) + sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = ONE + lambda*(X1 - X0);
  w2 = ONE + lambda*(X2 - X0);
  w3 = ONE + lambda*(X3 - X0);
  w4 = ONE + lambda*(X4 - X0);
  // Determine the interpolated state:
  return (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Ww?? -- Cell node primitive solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwNW(const int ii, const int jj) {
  return Ww(ii,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwNE(const int ii, const int jj) {
  return Ww(ii+1,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwSE(const int ii, const int jj) {
  return Ww(ii+1,jj);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwSW(const int ii, const int jj) {
  return Ww(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Uw?? -- Cell node conserved solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwNW(const int ii, const int jj) {
  return Uw(ii,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwNE(const int ii, const int jj) {
  return Uw(ii+1,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwSE(const int ii, const int jj) {
  return Uw(ii+1,jj);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwSW(const int ii, const int jj) {
  return Uw(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::evaluate_limiters -- Set flag to        *
 *                                                 evaluate limiters. *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::freeze_limiters -- Set flag to freeze   *
 *                                               limiters.            *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
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
inline void NavierStokes2D_Quad_Block::SetPiecewiseLinearReconstructionStencil(const int &i, const int &j,
									       int i_index[], int j_index[],
									       int &n_pts){

  if (i <= ICl-Nghost || i >= ICu+Nghost || j <= JCl-Nghost || j >= JCu+Nghost) {
    n_pts = 0;
  } else if ((i == ICl-1) && (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
	       Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_INFLOW_SUBSONIC ||
	       Grid.BCtypeW[j] == BC_OUTFLOW_SUBSONIC ||
	       Grid.BCtypeW[j] == BC_FIXED_PRESSURE ||
	       Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
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
	       Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_INFLOW_SUBSONIC ||
	       Grid.BCtypeE[j] == BC_OUTFLOW_SUBSONIC ||
	       Grid.BCtypeE[j] == BC_FIXED_PRESSURE ||
	       Grid.BCtypeE[j] == BC_CHARACTERISTIC ||	       
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
	       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_INFLOW_SUBSONIC ||
	       Grid.BCtypeS[i] == BC_OUTFLOW_SUBSONIC ||
	       Grid.BCtypeS[i] == BC_FIXED_PRESSURE ||
	       Grid.BCtypeS[i] == BC_CHARACTERISTIC ||	       
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
	       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_INFLOW_SUBSONIC ||
	       Grid.BCtypeN[i] == BC_OUTFLOW_SUBSONIC ||
	       Grid.BCtypeN[i] == BC_FIXED_PRESSURE ||
	       Grid.BCtypeN[i] == BC_CHARACTERISTIC ||	       
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

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Input-output operators.               *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const NavierStokes2D_Quad_Block &SolnBlk) {
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << endl;
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << endl;
  out_file << SolnBlk.Nghost << endl;
  out_file << SolnBlk.Axisymmetric << endl;
  out_file << SolnBlk.Compressibility_Effect << endl;
  out_file << SolnBlk.Transition_Model << endl;
  out_file << SolnBlk.Variable_Prandtl << endl;
  out_file << SolnBlk.Flow_Type << endl;
  out_file << SolnBlk.Vwall << endl;
  out_file << SolnBlk.Twall << endl;
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return out_file;
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      out_file << SolnBlk.U[i][j];
      out_file << endl;
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    out_file << SolnBlk.WoW[j] << endl;
    out_file << SolnBlk.WoE[j] << endl;
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    out_file << SolnBlk.WoS[i] << endl;
    out_file << SolnBlk.WoN[i] << endl;
  }

  // Output the high-order variables
  out_file << SolnBlk.NumberOfHighOrderVariables << "\n";
  for (int k = 1; k <= SolnBlk.NumberOfHighOrderVariables; ++k){
    out_file << SolnBlk.HighOrderVariable(k-1);
  }/* endfor */

  return out_file;
}

inline istream &operator >> (istream &in_file,
			     NavierStokes2D_Quad_Block &SolnBlk) {
  int i,j,k;
  int ni, il, iu, nj, jl, ju, ng, n_HO;
  Grid2D_Quad_Block_HO New_Grid;
  in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Compressibility_Effect;
  in_file >> SolnBlk.Transition_Model;
  in_file >> SolnBlk.Variable_Prandtl;
  in_file >> SolnBlk.Flow_Type;
  in_file.unsetf(ios::skipws);
  in_file >> SolnBlk.Vwall;
  in_file.setf(ios::skipws);
  in_file >> SolnBlk.Twall;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
    SolnBlk.deallocate(); return in_file;
  }
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) {
    if (SolnBlk.U != NULL) SolnBlk.deallocate();
    SolnBlk.allocate(ni-2*ng,nj-2*ng,ng);
  }

  // Copy the temporary mesh into the grid of the current solution block
  SolnBlk.Grid = New_Grid;

  // Read the solution & Initialize some data structures
  for (j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      in_file >> SolnBlk.U[i][j];
      SolnBlk.W[i][j] = SolnBlk.U[i][j].W();
      for (k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; k++) {
	SolnBlk.dUdt[i][j][k].Vacuum();
      }
      SolnBlk.dWdx[i][j].Vacuum();
      SolnBlk.dWdy[i][j].Vacuum();
      SolnBlk.phi[i][j].Vacuum();
      SolnBlk.Uo[i][j].Vacuum();
      SolnBlk.dt[i][j] = ZERO;
    }
  }
  if (SolnBlk.face_grad_arrays_allocated) {
    for (int i = 0; i < SolnBlk.NCi; i++) { 
      for (int j = 0; j < SolnBlk.NCj; j++) {
        SolnBlk.dWdx_faceN[i][j].Vacuum(); SolnBlk.dWdy_faceN[i][j].Vacuum();
        SolnBlk.dWdx_faceE[i][j].Vacuum(); SolnBlk.dWdy_faceE[i][j].Vacuum();
      }
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    in_file >> SolnBlk.WoW[j];
    in_file >> SolnBlk.WoE[j];
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    in_file >> SolnBlk.WoS[i];
    in_file >> SolnBlk.WoN[i];
  }

  // Allocate memory for the high-order boundary conditions (the values don't need to be read)
  SolnBlk.allocate_HighOrder_BoundaryConditions();

  in_file.setf(ios::skipws);

  // Read the high-order variables
  in_file >> n_HO;
  SolnBlk.allocate_HighOrder_Array(n_HO);

  for (int k = 1; k <= SolnBlk.NumberOfHighOrderVariables; ++k){
    // Read the current variable
    in_file >> SolnBlk.HighOrderVariable(k-1);
    // Associate this variable to the current grid
    SolnBlk.HighOrderVariable(k-1).AssociateGeometry(SolnBlk.Grid);
  }/* endfor */

  return in_file;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::NumVar -- Returns number of state       *
 *                                      variables.                    *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::NumVar(void) {
  return NUM_VAR_NAVIERSTOKES2D;
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
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
											const double &DeltaXToCentroid,
											const double &DeltaYToCentroid) const{
  return W[ii][jj] + (phi[ii][jj]^(dWdx[ii][jj]*DeltaXToCentroid + dWdy[ii][jj]*DeltaYToCentroid) );
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
inline double NavierStokes2D_Quad_Block::PiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
									 const double &DeltaXToCentroid,
									 const double &DeltaYToCentroid,
									 const unsigned int &parameter) const{
  return W[ii][jj][parameter] + (phi[ii][jj][parameter]*(dWdx[ii][jj][parameter]*DeltaXToCentroid + 
							 dWdy[ii][jj][parameter]*DeltaYToCentroid) );
}

/*!
 * Compute the solution provided by the piecewise linear reconstruction
 * at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline NavierStokes2D_pState 
NavierStokes2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
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
inline double NavierStokes2D_Quad_Block::PiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
									   const Vector2D &CalculationPoint,
									   const unsigned int &parameter) const{
  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	  
  return PiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y,parameter);
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
inline NavierStokes2D_pState 
NavierStokes2D_Quad_Block::UnlimitedPiecewiseLinearSolutionForDelta(const int &ii, const int &jj,
								    const double &DeltaXToCentroid,
								    const double &DeltaYToCentroid) const{

  return W[ii][jj] + dWdx[ii][jj]*DeltaXToCentroid + dWdy[ii][jj]*DeltaYToCentroid;
}

/*!
 * Compute the solution provided by the UNLIMITED (i.e. no limiter!) piecewise 
 * linear reconstruction at a particular location.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 */
inline NavierStokes2D_pState 
NavierStokes2D_Quad_Block::UnlimitedPiecewiseLinearSolutionAtLocation(const int &ii, const int &jj,
								      const Vector2D &CalculationPoint) const{

  // calculate the distance between the point of interest and the centroid of the cell
  Vector2D dX(CalculationPoint - Grid.Cell[ii][jj].Xc);	
  return UnlimitedPiecewiseLinearSolutionForDelta(ii,jj,dX.x,dX.y);
}

/*!
 * Compute the solution entropy at a particular location 
 * based on the piecewise linear reconstruction.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint the Cartesian coordinates of the point of interest
 * \param parameter only used for compatibility with other frameworks, otherwise useless
 */
inline double NavierStokes2D_Quad_Block::SolutionEntropyAtLocation(const int &ii, const int &jj,
								   const Vector2D &CalculationPoint,
								   const unsigned int &parameter) const{
  
  return PiecewiseLinearSolutionAtLocation(ii,jj,CalculationPoint).s();
}

/*!
 * Compute the solution pressure at a particular location 
 * based on the piecewise linear reconstruction of cell (ii,jj).
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param X_Coord the x-coordinate of the point of interest
 * \param Y_Coord the y-coordinate of the point of interest
 */
inline double NavierStokes2D_Quad_Block::SolutionPressureAtCoordinates(const int &ii, const int &jj,
								       const double &X_Coord, const double &Y_Coord) const{

  return W[ii][jj][4] + (phi[ii][jj][4]*(dWdx[ii][jj][4]*(X_Coord - Grid.XCellCentroid(ii,jj)) + 
					 dWdy[ii][jj][4]*(Y_Coord - Grid.YCellCentroid(ii,jj)) ));
}

/*!
 * Compute the wall shear stress at a particular location
 * based on the high-order piecewise reconstruction of cell (ii,jj).
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param CalculationPoint position vector for the location of interest
 * \param normal_dir normal vector at the location of interest
 */
inline double NavierStokes2D_Quad_Block::WallShearStress_HighOrder(const int &ii, const int &jj,
								   const Vector2D &CalculationPoint,
								   const Vector2D &normal_dir) const {
  // Solution and gradient in the normal direction at the location of interest
  NavierStokes2D_pState W_P;
  double dudn_P, dvdn_P;

  // == determine solution at the given location
  W_P = HO_Ptr[0].SolutionStateAtLocation(ii,jj,CalculationPoint);
  // == determine normal gradient of x-velocity component 
  dudn_P = HO_Ptr[0].NormalGradientAtLocation(ii,jj,
					      CalculationPoint, normal_dir,
					      2); // < x-velocity component
  // == determine normal gradient of y-velocity component 
  dvdn_P = HO_Ptr[0].NormalGradientAtLocation(ii,jj,
					      CalculationPoint, normal_dir,
					      3); // < y-velocity component

  // return the value of the wall shear stress (i.e a scalar)
  return W_P.WallShearStress(dudn_P, dvdn_P, normal_dir);
}

/*!
 * Allocate the high-order array.
 */
inline void NavierStokes2D_Quad_Block::allocate_HighOrder_Array(const int & NumberOfReconstructions){

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
 * Copy the high-order objects from the SolnBlk.
 */
inline void NavierStokes2D_Quad_Block::copy_HighOrder_Objects(const NavierStokes2D_Quad_Block &SolnBlk){

  int i, j, k;

  if (SolnBlk.U != NULL){

    /* Set the same number of high-order objects
       as that of the rhs block. */
    allocate_HighOrder_Array(SolnBlk.NumberOfHighOrderVariables);
    
    // Copy the high-order objects
    for (k = 1; k <= NumberOfHighOrderVariables; ++k){
      // set geometry pointer for each high-order object to the current grid
      HighOrderVariable(k-1).SetGeometryPointer(Grid);
      // copy high-order content from the the SolnBlk
      HighOrderVariable(k-1) = SolnBlk.HighOrderVariable(k-1);
    }/* endfor */

    // Copy the high-order boundary conditions.

    // === North BCs
    if (SolnBlk.HO_WoN != NULL){

      if (HO_WoN != NULL){
	// deallocate memory
	delete [] HO_WoN; HO_WoN = NULL;
      }

      // allocate new memory based on the number of the current grid
      HO_WoN = new BC_Type[NCi];

      for (i=0; i<NCi; ++i){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_NorthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_North(i,JCu),
					    Grid.BCtypeN[i]);

	// Copy North high-order BCs
	HO_WoN[i] = SolnBlk.HO_WoN[i];
      }
      
    } else if ( HO_WoN != NULL){
      // deallocate memory
      delete [] HO_WoN; HO_WoN = NULL;
    }


    // === South BCs
    if (SolnBlk.HO_WoS != NULL){

      if (HO_WoS != NULL){
	// deallocate memory
	delete [] HO_WoS; HO_WoS = NULL;
      }

      // allocate new memory    
      HO_WoS = new BC_Type[NCi];

      for (i=0; i<NCi; ++i){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_SouthCell(i).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_South(i,JCl),
					    Grid.BCtypeS[i]);
	
	// Copy South high-order BCs
	HO_WoS[i] = SolnBlk.HO_WoS[i];
      }

    } else if (HO_WoS != NULL){
      // deallocate memory
      delete [] HO_WoS; HO_WoS = NULL;
    }


    // === East BCs
    if (SolnBlk.HO_WoE != NULL){

      if (HO_WoE != NULL){
	// deallocate memory
	delete [] HO_WoE; HO_WoE = NULL;
      }

      // allocate new memory    
      HO_WoE = new BC_Type[NCj];

      for (j=0; j<NCj; ++j){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_EastCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_East(ICu,j),
					   Grid.BCtypeE[j]);

	// Copy East high-order BCs
	HO_WoE[j] = SolnBlk.HO_WoE[j];
      }
    } else if (HO_WoE != NULL){
      // deallocate memory
      delete [] HO_WoE; HO_WoE = NULL;
    }

    // === West BCs
    if (SolnBlk.HO_WoW != NULL){

      if (HO_WoW != NULL){
	// deallocate memory
	delete [] HO_WoW; HO_WoW = NULL;
      }

      // allocate new memory    
      HO_WoW = new BC_Type[NCj];
      
      for (j=0; j<NCj; ++j){
	// allocate BC memory for each constrained Gauss quadrature point
	BC_WestCell(j).InitializeCauchyBCs(Grid.NumOfConstrainedGaussQuadPoints_West(ICl,j),
					   Grid.BCtypeW[j]);
	
	// Copy West high-order BCs
	HO_WoW[j] = SolnBlk.HO_WoW[j];
      }

    } else if (HO_WoW != NULL){
      // deallocate memory
      delete [] HO_WoW; HO_WoW = NULL;
    }
    
  } // endif
  
}


/**********************************************************************
 *                                                                    *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                     *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer -- Loads send message    *
 *                                              buffer.               *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer(double *buffer,
						     int &buffer_count,
						     const int buffer_size,
						     const int i_min,
						     const int i_max,
						     const int i_inc,
						     const int j_min,
						     const int j_max,
						     const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = U[i][j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer_F2C --                   *
 *                     Loads send message buffer for fine to coarse   *
 *                     block message passing.                         *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
							 int &buffer_count,
							 const int buffer_size,
							 const int i_min,
							 const int i_max,
							 const int i_inc,
							 const int j_min,
							 const int j_max,
							 const int j_inc) {

  for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
				Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
				Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
				Grid.Cell[i+1][j+1].A*U[i+1][j+1][k])/
	                       (Grid.Cell[i  ][j  ].A+
				Grid.Cell[i+1][j  ].A+
				Grid.Cell[i  ][j+1].A+
				Grid.Cell[i+1][j+1].A);
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive  *
 *                                                   message buffer.  *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
							  int &buffer_count,
							  const int buffer_size,
							  const int i_min,
							  const int i_max,
							  const int i_inc,
							  const int j_min,
							  const int j_max,
							  const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_F2C --              *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing.                  *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
							      int &buffer_count,
							      const int buffer_size,
							      const int i_min,
							      const int i_max,
							      const int i_inc,
							      const int j_min,
							      const int j_max,
							      const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;    
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_C2F --              *
 *                     Unloads receive message buffer for coarse to   *
 *                     fine block message passing.                    *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
							      int &buffer_count,
							      const int buffer_size,
							      const int i_min,
							      const int i_max,
							      const int i_inc,
							      const int j_min,
							      const int j_max,
							      const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::SubcellReconstruction --                *
 *                     Performs the subcell linear reconstruction of  *
 *                     solution state within a given cell (i,j) of    *
 *                     the computational mesh for the specified       *
 *                     quadrilateral solution block.                  *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::SubcellReconstruction(const int i,
							     const int j,
							     const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  NavierStokes2D_pState DU, DUDx_ave, DUDy_ave;

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.

  if (i == ICl-2 || i == ICu+2 || j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if (i == ICl-1 && Grid.BCtypeW[j] != BC_NONE) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
               Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
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
      }
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
       }
    }
  } else if (i == ICu+1 && Grid.BCtypeE[j] != BC_NONE) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
               Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
  } else if (j == JCl-1 && Grid.BCtypeS[i] != BC_NONE) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
               Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
  } else if (j == JCu+1 && Grid.BCtypeN[i] != BC_NONE) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
               Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
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
  }
  
  if (n_pts > 0) {
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
  
    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = Grid.Cell[i_index[n2]][j_index[n2]].Xc - Grid.Cell[i][j].Xc;
      DU = W[i_index[n2]][j_index[n2]] - W[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    }
  					    
    DUDx_ave /= double(n_pts);
    DUDy_ave /= double(n_pts);
    DxDx_ave /= double(n_pts);
    DxDy_ave /= double(n_pts);
    DyDy_ave /= double(n_pts);
    dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
    // Calculate slope limiters. 
    if (!Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR_NAVIERSTOKES2D; n++) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < n_pts; n2++) {
	  u0Min = min(u0Min,W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,W[i_index[n2]][j_index[n2]][n]);
	}

	dX = Grid.xfaceE(i,j) - Grid.Cell[i][j].Xc;
	uQuad[0] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceW(i,j) - Grid.Cell[i][j].Xc;
	uQuad[1] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceN(i,j) - Grid.Cell[i][j].Xc;
	uQuad[2] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceS(i,j) - Grid.Cell[i][j].Xc;
	uQuad[3] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
    
	switch(Limiter) {
	case LIMITER_ONE :
	  phi_n = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_n = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_n = Limiter_BarthJespersen(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_n = Limiter_Venkatakrishnan(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANLEER :
	  phi_n = Limiter_VanLeer(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANALBADA :
	  phi_n = Limiter_VanAlbada(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	default:
	  phi_n = Limiter_BarthJespersen(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	}

	phi[i][j][n] = phi_n;

      }
    }
  } else {
    dWdx[i][j].Vacuum();
    dWdy[i][j].Vacuum();
    phi[i][j].Vacuum();
  }

}

/*******************************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer_Flux_F2C --                       *
 *                     Loads send message buffer for fine to coarse block      *
 *                     message passing of conservative solution fluxes.        *
 *******************************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
							      int &buffer_count,
							      const int buffer_size,
							      const int i_min,
							      const int i_max,
							      const int i_inc,
							      const int j_min,
							      const int j_max,
							      const int j_inc) {

  if (j_min == j_max && j_min == JCl) {
    for (int i = i_min;  ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxS[i  ][k] + FluxS[i+1][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxN[i  ][k] + FluxN[i+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxW[j][k] + FluxW[j+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxE[j][k] + FluxE[j+1][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C --         *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing of conservative   *
 *                     solution fluxes.                               *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
								   int &buffer_count,
								   const int buffer_size,
								   const int i_min,
								   const int i_max,
								   const int i_inc,
								   const int j_min,
								   const int j_max,
								   const int j_inc) {

  if (j_min == j_max && j_min == JCl) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    }
  }
  return 0;
}

/*!
 * Return the upwind inviscid flux in the normal direction 
 * based on the left and right interface states for 
 * a variety of flux functions.
 *
 * \param Flux_Function index to specify the requested flux function
 * \param Wl left interface state
 * \param Wr right interface state
 * \param normal_dir vector to define the normal direction
 */
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::RiemannFlux_n(const int & Flux_Function,
								      const NavierStokes2D_pState &Wl,
								      const NavierStokes2D_pState &Wr,
								      const Vector2D &normal_dir) const{

  switch(Flux_Function) {
  case FLUX_FUNCTION_GODUNOV :
    return FluxGodunov_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_ROE :
    return FluxRoe_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_RUSANOV :
    return FluxRusanov_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_HLLE :
    return FluxHLLE_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_HLLL :
    return FluxHLLL_n(Wl,Wr, normal_dir);
  case FLUX_FUNCTION_HLLC :
    return FluxHLLC_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_VANLEER :
    return FluxVanLeer_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_AUSM :
    return FluxAUSM_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_AUSMplus :
    return FluxAUSMplus_n(Wl, Wr, normal_dir);
  case FLUX_FUNCTION_GODUNOV_WRS :
    return FluxGodunov_Wrs_n(Wl,Wr, normal_dir);
  default:
    return FluxRoe_n(Wl, Wr, normal_dir);
  } /* endswitch */

}

/*!
 * Check that the primitive state is physical.
 *
 * \param W checked solution state
 * \param iCell i-index of the cell which the solution state belongs to 
 * \param iCell j-index of the cell which the solution state belongs to 
 * \param Ref reference string. Used in the output error message
 */
inline void NavierStokes2D_Quad_Block::Validate_Primitive_SolnState(NavierStokes2D_pState & W,
								    const int &iCell,
								    const int &jCell,
								    const std::string &Ref,
								    const int &IndexHO) const{
  
  int Param_Index;
  std::ostringstream error_msg;

  // Check if negative density or pressure occur
  if ( W.Unphysical_Properties() ) {

    if (CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE){

      if (CENO_Execution_Mode::CENO_VERBOSE){ 
	// output a brief error message
	std::cout << "\n " << CFFC_Name() 
		  << " NavierStokes2D ERROR: Unphysical primitive variable at the "
		  << Ref
		  << " interface of (" << iCell << "," << jCell << ") cell: "
		  << "\n W = " << W << "\n"; 
      }

      // try using the Piecewise Constant (PWC) solution instead
      W = CellSolution(iCell,jCell);

    } else {

      // throw a runtime error with an error message
      error_msg << "\n " << CFFC_Name() 
		<< " NavierStokes2D ERROR: Unphysical primitive variable at the "
		<< Ref
		<< " interface of (" << iCell << "," << jCell << ") cell: "
		<< "\n W = " << W; 
      
      error_msg << "\n Cell Centroid: "
		<< Grid.Cell[iCell][jCell].Xc;

      error_msg << "\n High-order reconstruction data: "
		<< "\n Derivatives: \n"
		<< HighOrderVariable(IndexHO).CellTaylorDeriv(iCell,jCell)
		<< "\n Limiter: \n"
		<< HighOrderVariable(IndexHO).CellTaylorDeriv(iCell,jCell).Limiter();
	      
      for (Param_Index = 1; Param_Index <= 4; ++Param_Index){
	error_msg << "\n  Variable=" << Param_Index
		  << ", SI=" << HighOrderVariable(IndexHO).CellSmoothnessIndicatorValue(iCell,jCell,Param_Index)
		  << ", Flagged=" <<  HighOrderVariable(IndexHO).CellInadequateFitValue(iCell,jCell,Param_Index);
      }
      error_msg << "\n";

      throw runtime_error(error_msg.str());

    } // endif (FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE)

  } // endif

}

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Single Block External Subroutines.    *
 **********************************************************************/

extern void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator,
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk1,
		                NavierStokes2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Fine,
				  NavierStokes2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Coarse,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_SW,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_SE,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_NW,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_NE);

extern void ICs(NavierStokes2D_Quad_Block &SolnBlk,
		NavierStokes2D_Input_Parameters &IP,
                NavierStokes2D_pState *Wo);

extern void BCs(NavierStokes2D_Quad_Block &SolnBlk,
		NavierStokes2D_Input_Parameters &IP);

extern double CFL(NavierStokes2D_Quad_Block &SolnBlk,
                  NavierStokes2D_Input_Parameters &IP);

extern void Set_Global_TimeStep(NavierStokes2D_Quad_Block &SolnBlk,
                                const double &Dt_min);

extern double L1_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);
 
extern double L1_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int residual_var);

extern double L2_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int residual_var);

extern double Max_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk, int residual_var);

extern void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(NavierStokes2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction(NavierStokes2D_Quad_Block &SolnBlk,
				  const int & Reconstruction_Type,
				  const int & Limiter);

extern void Residual_Smoothing(NavierStokes2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon,
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  NavierStokes2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          NavierStokes2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(NavierStokes2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
								NavierStokes2D_Input_Parameters &IP,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(NavierStokes2D_Quad_Block &SolnBlk,
				    NavierStokes2D_Input_Parameters &IP);

extern int dUdt_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    NavierStokes2D_Input_Parameters &IP);

extern int Update_Solution_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
					       NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Multiple Block External Subroutines.  *
 **********************************************************************/

extern NavierStokes2D_Quad_Block* Allocate(NavierStokes2D_Quad_Block *Soln_ptr,
					   NavierStokes2D_Input_Parameters &Input_Parameters);

extern NavierStokes2D_Quad_Block* Deallocate(NavierStokes2D_Quad_Block *Soln_ptr,
					     NavierStokes2D_Input_Parameters &Input_Parameters);

extern NavierStokes2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block_HO **InitMeshBlk,
							      NavierStokes2D_Quad_Block *Soln_ptr,
							      NavierStokes2D_Input_Parameters &Input_Parameters,
							      QuadTreeBlock_DataStructure &QuadTree,
							      AdaptiveBlockResourceList &GlobalSolnBlockList,
							      AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(NavierStokes2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                NavierStokes2D_Input_Parameters &Input_Parameters);

extern void Linear_Reconstruction(NavierStokes2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  NavierStokes2D_Input_Parameters &Input_Parameters);

extern void BCs(NavierStokes2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		NavierStokes2D_Input_Parameters &Input_Parameters);

extern double CFL(NavierStokes2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
		  NavierStokes2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                const double &Dt_min);

extern double L1_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

void L1_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List, double *array);

void L2_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List, double *array);

void Max_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, double *array);

extern void Evaluate_Limiters(NavierStokes2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(NavierStokes2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(NavierStokes2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                NavierStokes2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Residual_Evaluation(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    NavierStokes2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    NavierStokes2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               NavierStokes2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
				    NavierStokes2D_Input_Parameters &Input_Parameters,
				    QuadTreeBlock_DataStructure &QuadTree,
				    AdaptiveBlockResourceList &GlobalSolnBlockList,
				    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Initial_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					    NavierStokes2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Uniform_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					    NavierStokes2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Boundary_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					     NavierStokes2D_Input_Parameters &Input_Parameters,
					     QuadTreeBlock_DataStructure &QuadTree,
					     AdaptiveBlockResourceList &GlobalSolnBlockList,
					     AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Flat_Plate_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					       NavierStokes2D_Input_Parameters &Input_Parameters,
					       QuadTreeBlock_DataStructure &QuadTree,
					       AdaptiveBlockResourceList &GlobalSolnBlockList,
					       AdaptiveBlock2D_List &LocalSolnBlockList);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- IO Single Block External Subroutines. *
 **********************************************************************/

extern void Write_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Output_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			   NavierStokes2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		 	         NavierStokes2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nozzleless_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				      const int Number_of_Time_Steps,
				      const double &Time,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File,
				      const double &po,
				      const double &rhoo);

extern void Output_Nodes_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Gradients_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Quasi3D_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				   NavierStokes2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

extern void Output_Ringleb_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   double &area);

extern void Output_Viscous_Channel_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					   NavierStokes2D_Input_Parameters &IP,
					   const int Block_Number,
					   const int Output_Title,
					   ostream &Out_File,
					   double &l1_norm,
					   double &l2_norm,
					   double &max_norm);

extern void Output_Viscous_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					NavierStokes2D_Input_Parameters &IP,
					const int Block_Number,
					const int Output_Title,
					ostream &Out_File,
					double &l1_norm,
					double &l2_norm,
					double &max_norm);

extern void Output_Turbulent_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					  const int Block_Number,
					  const int Output_Title,
					  const int Output_Data,
					  ostream &Out_File,
					  const double &Re,
					  const double &Pipe_Radius,
					  const int &variable_flag);

extern void Output_Subsonic_Hot_Jet_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					    const int Block_Number,
					    const int Output_Title,
					    const int Output_Data,
					    ostream &Out_File,
					    const int &variable_flag);

extern void Output_Supersonic_Hot_Jet_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					      const int Block_Number,
					      const int Output_Title,
					      const int Output_Data,
					      ostream &Out_File,
					      const int &variable_flag);

extern void Output_Flat_Plate_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				      const int Block_Number,
				      const int Output_Title_Soln,
				      ostream &Out_File_Soln,
				      const int Output_Title_Skin,
				      ostream &Out_File_Skin,
				      const NavierStokes2D_pState &Winf,
				      const double &plate_length,
				      double &l1_norm,
				      double &l2_norm,
				      double &max_norm,
				      double &area,
				      int &numberofactivecells,
				      double &l1_norm_cf,
				      double &l2_norm_cf,
				      double &max_norm_cf,
				      double &area_cf,
				      int &numberofactivecells_cf);

extern void Output_Driven_Cavity_Flow_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					      const int Block_Number,
					      const int Output_Title,
					      ostream &Out_File_u,
					      ostream &Out_File_v,
					      const double &Re,
					      const Vector2D &Vwall,
					      const double &length);

extern void Output_Backward_Facing_Step_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
						const int Block_Number,
						const int Output_Title,
						ostream &Out_File,
						const double &step_height,
						const double &top_wall_deflection);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- IO Multiple Block External Subroutines.*
 **********************************************************************/

extern int Read_Restart_Solution(NavierStokes2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 NavierStokes2D_Input_Parameters &IP,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(NavierStokes2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  NavierStokes2D_Input_Parameters &IP,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          NavierStokes2D_Input_Parameters &IP,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                NavierStokes2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                NavierStokes2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    NavierStokes2D_Input_Parameters &IP,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Quasi3D_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  NavierStokes2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  NavierStokes2D_Input_Parameters &IP);

extern int Output_Viscous_Channel_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &Soln_Block_List,
					  NavierStokes2D_Input_Parameters &IP);

extern int Output_Viscous_Pipe_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       NavierStokes2D_Input_Parameters &IP);

extern int Output_Turbulent_Pipe_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					 AdaptiveBlock2D_List &Soln_Block_List,
					 NavierStokes2D_Input_Parameters &IP);

extern int Output_Flat_Plate_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     NavierStokes2D_Input_Parameters &IP);

extern int Output_Supersonic_Hot_Jet_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					     AdaptiveBlock2D_List &Soln_Block_List,
					     NavierStokes2D_Input_Parameters &IP);

extern int Output_Subsonic_Hot_Jet_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					   AdaptiveBlock2D_List &Soln_Block_List,
					   NavierStokes2D_Input_Parameters &IP);

extern int Output_Driven_Cavity_Flow_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					     AdaptiveBlock2D_List &Soln_Block_List,
					     NavierStokes2D_Input_Parameters &IP);

extern int Output_Backward_Facing_Step_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					       AdaptiveBlock2D_List &Soln_Block_List,
					       NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Grid Multiple Block External          *
 *                              Subroutines.                          *
 **********************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					    NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						      NavierStokes2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
							    NavierStokes2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						 NavierStokes2D_Input_Parameters &IP);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          NavierStokes2D_Input_Parameters &IP);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                NavierStokes2D_Input_Parameters &IP);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Single Block External      *
 *                              Subroutines.                          *
 **********************************************************************/

extern int Turbulence_ICs(NavierStokes2D_Quad_Block &SolnBlk,
			  NavierStokes2D_Input_Parameters &IP,
			  NavierStokes2D_pState *Wo);

extern int Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
			  NavierStokes2D_Input_Parameters &IP);

extern int Apply_Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
				NavierStokes2D_Input_Parameters &IP,
				const int &Turbulent_BCtype,
				const int &i, const int &j);

extern int Turbulence_Zero_Residual(NavierStokes2D_Quad_Block &SolnBlk,
				    const int i_stage,
				    NavierStokes2D_Input_Parameters &IP);

extern int Zero_Turbulence_Residuals(NavierStokes2D_Quad_Block &SolnBlk,
				     NavierStokes2D_Input_Parameters &IP,
				     const int &Turbulent_BCtype,
				     const int &i, const int &j,
				     const int &k_residual);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Multiple Block External    *
 *                              Subroutines.                          *
 **********************************************************************/

extern int Turbulence_BCs(NavierStokes2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  NavierStokes2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Solvers.                              *
 **********************************************************************/

extern int NavierStokes2DQuadSolver(char *Input_File_Name_ptr,
				    int batch_flag);

/*
 * Include specializations CFFC header files.
 * Must be included at the end of the file!!!
 */
#include "NavierStokes2DAccuracyAssessment.h" /* Include 2D accuracy assessment header 
						 file with specializations for Navier-Stokes solution. */
#include "NavierStokes2DHighOrder.h"   /* Include 2D high-order header file with specializations for Navier-Stokes solution. */

#endif // _NAVIERSTOKES2D_QUAD_INCLUDED
